// Internal-only binary stage capture. Each group owns its stream; no shared
// append log, sequence copy, or production checkpoint is introduced.
struct DiagnosticSink {
    std::filesystem::path root;
    bool candidates_only{false};
    bool interval_trace{false};
    explicit DiagnosticSink(const std::string& path,bool only=false) : root(path),candidates_only(only) {
        const auto* format=std::getenv("RAMAG_COVERAGE_TRACE_FORMAT");
        if(format&&std::string_view(format)!="intervals-v1")
            throw std::runtime_error("unknown internal coverage trace format");
        interval_trace=format!=nullptr;
        if (!std::filesystem::create_directory(root))
            throw std::runtime_error("diagnostic directory must be new");
    }
    static void word(std::ostream& out,uint64_t value) {
        for(unsigned i=0;i<8;++i)out.put(static_cast<char>((value>>(8*i))&255));
    }
    static uint64_t word(std::istream& in) {
        uint64_t value=0;
        for(unsigned i=0;i<8;++i){const auto c=in.get();if(c==EOF)throw std::runtime_error("truncated diagnostic stream");value|=uint64_t(static_cast<unsigned char>(c))<<(8*i);}
        return value;
    }
    static std::ofstream output(const std::filesystem::path& p) {
        std::ofstream out(p,std::ios::binary|std::ios::out);
        out.exceptions(std::ios::badbit|std::ios::failbit);return out;
    }
    void seeds(std::span<const Seed> items) const {
        if(candidates_only||interval_trace)return; // input group trace already contains every seed
        auto out=output(root/"seeds.bin");out.write("RGS00001",8);word(out,items.size());
        for(const auto& s:items){word(out,s.reference_id);word(out,s.reference_begin);word(out,s.query_id);word(out,s.query_begin);word(out,s.length);word(out,s.strand==FORWARD?0:1);}
        out.flush();out.close();
    }
    void anchors(std::string_view stage,const AnchorVec& items) const {
        if(stage!="candidates") {
            auto out=output(root/(std::string(stage)+".bin"));out.write("RGF00001",8);word(out,items.size());
            for(size_t id=0;id<items.size();++id){word(out,id);word(out,items[id].ref_selected);word(out,items[id].qry_selected);}
            out.flush();out.close();
            return;
        }
        const auto path=root/(std::string(stage)+".bin"),partial=root/(std::string(stage)+".bin.partial");
        auto out=output(partial);out.write("RGA00002",8);word(out,items.size());
        for(size_t id=0;id<items.size();++id){const auto& a=items[id];word(out,id);word(out,a.reference_id);word(out,a.reference_begin);word(out,a.ref_len);word(out,a.query_id);word(out,a.query_begin);word(out,a.qry_len);word(out,a.strand==FORWARD?0:1);word(out,a.alignment_length);word(out,a.aligned_base);word(out,a.ref_selected);word(out,a.qry_selected);word(out,a.cigar.size());for(auto op:a.cigar)word(out,op);}
        out.write("RGAEND02",8);word(out,items.size());out.flush();out.close();
        auto sources=output(root/"candidates.sources.tsv.partial");sources<<"parent\tsource\n";
        for(size_t id=0;id<items.size();++id)sources<<id<<'\t'<<static_cast<unsigned>(items[id].candidate_source)<<'\n';
        sources<<"END\t"<<items.size()<<'\n';sources.flush();sources.close();
        std::filesystem::rename(root/"candidates.sources.tsv.partial",root/"candidates.sources.tsv");
        std::filesystem::rename(partial,path);
    }
    void output_sources(std::span<const PairwiseAlignment> records) const {
        auto out=output(root/"output.sources.tsv.partial");
        out<<"record\treference\trbegin\trend\tquery\tqbegin\tqend\treverse\tsource\tselected\n";
        for(size_t i=0;i<records.size();++i){const auto& x=records[i];const auto& r=x.record;
            out<<i<<'\t'<<r.reference_id<<'\t'<<r.reference_begin<<'\t'<<r.reference_end<<'\t'
               <<r.query_id<<'\t'<<r.query_begin<<'\t'<<r.query_end<<'\t'<<(r.strand==REVERSE)<<'\t'
               <<static_cast<unsigned>(x.candidate_source)<<'\t'<<(x.reference_selected&&x.query_selected)<<'\n';}
        out<<"END\t"<<records.size()<<'\n';out.flush();out.close();
        std::filesystem::rename(root/"output.sources.tsv.partial",root/"output.sources.tsv");
    }
};
thread_local std::ostream* diagnostic_group_stream=nullptr;
// Fixed-size records, no per-seed text formatting or all-stage resident copy.
// Fields: stage, parent, source, reference ID/begin/end, query ID/begin/end.
thread_local std::ostream* diagnostic_interval_stream=nullptr;
void TraceInterval(uint64_t stage,uint64_t parent,uint64_t source,uint64_t rid,
                   uint64_t rb,uint64_t re,uint64_t qid,uint64_t qb,uint64_t qe) {
    if(!diagnostic_interval_stream)return;
    for(auto value:{stage,parent,source,rid,rb,re,qid,qb,qe})
        DiagnosticSink::word(*diagnostic_interval_stream,value);
}
uint64_t TraceStage(std::string_view stage) {
    constexpr std::string_view names[]={"input","merged","cluster","best-chain",
        "supported-chain","pre-cleanup","post-cleanup","extended","linked","supplemental"};
    for(uint64_t i=0;i<std::size(names);++i)if(stage==names[i])return i;
    throw std::runtime_error("unregistered coverage interval stage");
}
struct DiagnosticGroupScope {
    std::ofstream stream;
    std::ofstream intervals;
    std::filesystem::path interval_path;
    std::ostream* previous=diagnostic_group_stream;
    std::ostream* previous_intervals=diagnostic_interval_stream;
    DiagnosticGroupScope(const DiagnosticSink* sink,size_t group,std::string_view phase) {
        if(sink&&(!sink->candidates_only||phase=="selection")){
            if(sink->interval_trace&&phase!="selection") {
                interval_path=sink->root/(std::string(phase)+"-"+std::to_string(group)+".intervals");
                intervals=DiagnosticSink::output(interval_path.string()+".partial");
                intervals.write("RGI00001",8);diagnostic_interval_stream=&intervals;
            } else {stream=DiagnosticSink::output(sink->root/(std::string(phase)+"-"+std::to_string(group)+".tsv"));diagnostic_group_stream=&stream;}
        }
    }
    ~DiagnosticGroupScope(){diagnostic_group_stream=previous;diagnostic_interval_stream=previous_intervals;}
    void finish(){if(stream.is_open()){stream.flush();stream.close();}
        if(intervals.is_open()) {intervals.write("RGIEND01",8);intervals.flush();intervals.close();
            std::filesystem::rename(interval_path.string()+".partial",interval_path);}}
};
void TraceMatches(std::string_view stage,std::span<const Match> items,uint64_t parent=0) {
    if(diagnostic_interval_stream)for(const auto& s:items)
        TraceInterval(TraceStage(stage),parent,0,s.reference_id,s.reference_begin,s.reference_begin+s.length,
                      s.query_id,s.query_begin,s.query_begin+s.length);
    if(!diagnostic_group_stream)return;
    uint64_t id=0;
    for(const auto& s:items)*diagnostic_group_stream<<stage<<'\t'<<parent<<'\t'<<id++<<'\t'<<s.reference_id<<'\t'<<s.reference_begin<<'\t'<<s.query_id<<'\t'<<s.query_begin<<'\t'<<s.length<<'\t'<<(s.strand==FORWARD?'+':'-')<<'\n';
}
void TraceAnchors(std::string_view stage,const AnchorVec& items) {
    if(diagnostic_interval_stream)for(size_t id=0;id<items.size();++id) {
        const auto& a=items[id];uint64_t r=a.reference_begin,q=a.query_begin;
        if(a.strand==REVERSE)q+=a.qry_len;
        for(const auto unit:a.cigar) {
            const auto len=unit>>4,op=unit&15;
            if(op!=0&&op!=1&&op!=2&&op!=7&&op!=8)throw std::runtime_error("invalid diagnostic CIGAR");
            if(op!=1&&op!=2)TraceInterval(TraceStage(stage),id,a.candidate_source,a.reference_id,r,r+len,
                a.query_id,a.strand==REVERSE?q-len:q,a.strand==REVERSE?q:q+len);
            if(op!=1)r+=len;
            if(op!=2){if(a.strand==REVERSE)q-=len;else q+=len;}
        }
        if(r!=a.reference_begin+a.ref_len||q!=(a.strand==REVERSE?a.query_begin:a.query_begin+a.qry_len))
            throw std::runtime_error("diagnostic CIGAR consumption mismatch");
    }
    if(!diagnostic_group_stream)return;
    uint64_t id=0;
    for(const auto& a:items)*diagnostic_group_stream<<stage<<'\t'<<id++<<'\t'<<a.reference_id<<'\t'<<a.reference_begin<<'\t'<<a.ref_len<<'\t'<<a.query_id<<'\t'<<a.query_begin<<'\t'<<a.qry_len<<'\t'<<(a.strand==FORWARD?'+':'-')<<'\t'<<a.aligned_base<<'\t'<<a.alignment_length<<'\n';
}
// IDs are stable within the group/phase stream, even as a representative grows.
void TraceLinkDisposition(std::string_view event,size_t id,size_t parent,const Anchor& a) {
    if(!diagnostic_group_stream)return;
    *diagnostic_group_stream << "link-" << event << '\t' << id << '\t' << parent
        << '\t' << a.reference_id << '\t' << a.reference_begin << '\t' << a.ref_len
        << '\t' << a.query_id << '\t' << a.query_begin << '\t' << a.qry_len
        << '\t' << (a.strand==FORWARD?'+':'-') << '\t' << a.aligned_base
        << '\t' << a.alignment_length << '\n';
}

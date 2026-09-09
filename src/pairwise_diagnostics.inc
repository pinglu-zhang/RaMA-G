// Internal-only binary stage capture. Each group owns its stream; no shared
// append log, sequence copy, or production checkpoint is introduced.
struct DiagnosticSink {
    std::filesystem::path root;
    explicit DiagnosticSink(const std::string& path) : root(path) {
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
        auto out=output(root/(std::string(stage)+".bin"));out.write("RGA00001",8);word(out,items.size());
        for(size_t id=0;id<items.size();++id){const auto& a=items[id];word(out,id);word(out,a.reference_id);word(out,a.reference_begin);word(out,a.ref_len);word(out,a.query_id);word(out,a.query_begin);word(out,a.qry_len);word(out,a.strand==FORWARD?0:1);word(out,a.alignment_length);word(out,a.aligned_base);word(out,a.ref_selected);word(out,a.qry_selected);word(out,a.cigar.size());for(auto op:a.cigar)word(out,op);}
        out.flush();out.close();
    }
};
thread_local std::ostream* diagnostic_group_stream=nullptr;
struct DiagnosticGroupScope {
    std::ofstream stream;
    std::ostream* previous=diagnostic_group_stream;
    DiagnosticGroupScope(const DiagnosticSink* sink,size_t group,std::string_view phase) {
        if(sink){stream=DiagnosticSink::output(sink->root/(std::string(phase)+"-"+std::to_string(group)+".tsv"));diagnostic_group_stream=&stream;}
    }
    ~DiagnosticGroupScope(){diagnostic_group_stream=previous;}
    void finish(){if(stream.is_open()){stream.flush();stream.close();}}
};
void TraceMatches(std::string_view stage,std::span<const Match> items,uint64_t parent=0) {
    if(!diagnostic_group_stream)return;
    uint64_t id=0;
    for(const auto& s:items)*diagnostic_group_stream<<stage<<'\t'<<parent<<'\t'<<id++<<'\t'<<s.reference_id<<'\t'<<s.reference_begin<<'\t'<<s.query_id<<'\t'<<s.query_begin<<'\t'<<s.length<<'\t'<<(s.strand==FORWARD?'+':'-')<<'\n';
}
void TraceAnchors(std::string_view stage,const AnchorVec& items) {
    if(!diagnostic_group_stream)return;
    uint64_t id=0;
    for(const auto& a:items)*diagnostic_group_stream<<stage<<'\t'<<id++<<'\t'<<a.reference_id<<'\t'<<a.reference_begin<<'\t'<<a.ref_len<<'\t'<<a.query_id<<'\t'<<a.query_begin<<'\t'<<a.qry_len<<'\t'<<(a.strand==FORWARD?'+':'-')<<'\t'<<a.aligned_base<<'\t'<<a.alignment_length<<'\n';
}

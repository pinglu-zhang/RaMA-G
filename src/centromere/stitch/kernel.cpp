#include "centromere/stitch/kernel.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#if RAMAG_HAVE_CENTROMERE
#include "centrolign_stitch/alignment.hpp"
#endif

namespace ramag::detail {
namespace {
struct Scores {
  std::int64_t match, mismatch;
  std::array<std::int64_t,3> open, extend;
  std::size_t pieces;
};
Scores GetScores(CentromereStitchScoreModel model) {
  switch (model) {
    case CentromereStitchScoreModel::Centrolign3Piece:
      return {20,80,{60,800,2500},{30,5,1},3};
    case CentromereStitchScoreModel::Rama2PieceControl:
      return {0,3,{4,12,0},{2,1,0},2};
  }
  throw AlignmentError("unknown centromere stitch scoring model");
}
bool EqualBase(char a,char b) {
  return a==b && (a=='A'||a=='C'||a=='G'||a=='T');
}
void ValidateBase(char b) {
#if RAMAG_HAVE_CENTROMERE
  cg_budget::tick();
#endif
  if (b!='A'&&b!='C'&&b!='G'&&b!='T'&&b!='N')
    throw AlignmentError("centromere stitch requires uppercase A/C/G/T/N inputs");
}
void AddScore(std::int64_t& value,std::int64_t delta) {
  if ((delta>0&&value>std::numeric_limits<std::int64_t>::max()-delta)||
      (delta<0&&value<std::numeric_limits<std::int64_t>::min()-delta))
    throw AlignmentError("centromere stitch score overflow");
  value+=delta;
}
std::int64_t GapCost(std::uint64_t length,const Scores& p) {
  if (!length) return 0;
  auto cost=std::numeric_limits<std::int64_t>::max();
  for (std::size_t k=0;k<p.pieces;++k) {
    if (length>static_cast<std::uint64_t>((std::numeric_limits<std::int64_t>::max()-p.open[k])/p.extend[k]))
      throw AlignmentError("centromere stitch gap score overflow");
    cost=std::min(cost,p.open[k]+static_cast<std::int64_t>(length)*p.extend[k]);
  }
  return cost;
}
} // namespace

std::int64_t ScoreCentromereStitchCigar(std::string_view r,std::string_view q,
    std::span<const CigarOp> cigar,CentromereStitchScoreModel model) {
  const auto p=GetScores(model);
  std::size_t ri=0,qi=0;
  std::uint64_t gap_length=0;
  char gap_op=0;
  std::int64_t score=0;
  auto flush=[&](){AddScore(score,-GapCost(gap_length,p));gap_length=0;gap_op=0;};
  for (const auto& op:cigar) {
    if (!op.length) throw AlignmentError("zero-length centromere stitch CIGAR operation");
    if (op.operation!='='&&op.operation!='X'&&op.operation!='I'&&op.operation!='D')
      throw AlignmentError("invalid centromere stitch CIGAR operation");
    const bool rc=op.operation!='I',qc=op.operation!='D';
    if ((rc&&op.length>r.size()-ri)||(qc&&op.length>q.size()-qi))
      throw AlignmentError("centromere stitch CIGAR exceeds gap bounds");
    if (op.operation=='I'||op.operation=='D') {
      if (gap_op && gap_op!=op.operation) flush();
      gap_op=op.operation;
      if (op.length>std::numeric_limits<std::uint64_t>::max()-gap_length)
        throw AlignmentError("centromere stitch CIGAR length overflow");
      gap_length+=op.length;
      for (std::uint64_t i=0;i<op.length;++i) ValidateBase(rc?r[ri+i]:q[qi+i]);
    } else {
      flush();
      for (std::uint64_t i=0;i<op.length;++i) {
        ValidateBase(r[ri+i]);ValidateBase(q[qi+i]);
        const bool equal=EqualBase(r[ri+i],q[qi+i]);
        if ((op.operation=='=')!=equal)
          throw AlignmentError("centromere stitch EQ/X disagrees with FASTA or N mismatch rule");
        AddScore(score,equal?p.match:-p.mismatch);
      }
    }
    if(rc)ri+=op.length;
    if(qc)qi+=op.length;
  }
  flush();
  if(ri!=r.size()||qi!=q.size()) throw AlignmentError("centromere stitch CIGAR does not consume complete gap");
  return score;
}

#if RAMAG_HAVE_CENTROMERE
namespace {
namespace cv=centrolign_stitch_vendor;
char Encode(char base,bool query) {
  ValidateBase(base);
  switch(base) {case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:return query?5:4;}
}
cv::BaseGraph MakeLinearGraph(std::string_view sequence,bool query) {
  cv::BaseGraph g;
  for(std::size_t i=0;i<sequence.size();++i) {
    cg_budget::tick();g.add_node(Encode(sequence[i],query));if(i)g.add_edge(i-1,i);
  }
  return g;
}
// std::vector in the public result must keep its public allocator. Reserve its
// actual capacity under the same cap while it overlaps vendor graph/trace data.
struct OutputLease {
  cg_budget::Context& context;
  std::size_t bytes{};
  ~OutputLease(){context.release(bytes);}
  void reserve(std::vector<CigarOp>& v,std::size_t wanted) {
    if(wanted<=v.capacity())return;
    if(wanted>std::numeric_limits<std::size_t>::max()/sizeof(CigarOp))throw cg_budget::Limit{"gap_memory_budget"};
    const auto requested=wanted*sizeof(CigarOp);
    context.acquire(requested);
    try {v.reserve(wanted);}
    catch(const std::bad_alloc&){context.release(requested);throw cg_budget::Limit{"allocation_failure"};}
    catch(...){context.release(requested);throw;}
    const auto actual=v.capacity()*sizeof(CigarOp);
    if(actual>requested) {
      // Standard permits over-allocation; never silently account a smaller cap.
      try {context.acquire(actual-requested);}
      catch(...){context.release(requested);throw;}
    }
    context.release(bytes);bytes=actual;
  }
};
void Append(std::vector<CigarOp>& cigar,char op,OutputLease& lease) {
  if(!cigar.empty()&&cigar.back().operation==op){++cigar.back().length;return;}
  if(cigar.size()==cigar.capacity()) {
    const auto next=cigar.capacity()?cigar.capacity()*2:8;
    if(next<cigar.capacity())throw cg_budget::Limit{"gap_memory_budget"};
    lease.reserve(cigar,next);
  }
  cigar.push_back({op,1});
}
template<int Pieces>
StitchGapResult Run(std::string_view r,std::string_view q,const CentromereStitchOptions& options,
                   cg_budget::Context& context,bool dp,std::uint64_t cells) {
  const auto p=GetScores(options.score_model);
  cv::AlignmentParameters<Pieces> params;
  params.match=static_cast<std::uint32_t>(p.match);params.mismatch=static_cast<std::uint32_t>(p.mismatch);
  for(int i=0;i<Pieces;++i){params.gap_open[i]=static_cast<std::uint32_t>(p.open[i]);params.gap_extend[i]=static_cast<std::uint32_t>(p.extend[i]);}
  context.phase="stitch_gap_graphs";
  auto gr=MakeLinearGraph(r,false),gq=MakeLinearGraph(q,true);
  cg_budget::vector<std::uint64_t> sr{0},sq{0},er{r.size()-1},eq{q.size()-1};
  std::int64_t kernel_score=0;
  context.phase=dp?"stitch_gap_dp":"stitch_gap_wfa";
  auto alignment=dp?cv::po_poa<Pieces>(gr,gq,sr,sq,er,eq,params,&kernel_score)
                   :cv::wfa_po_poa<Pieces>(gr,gq,sr,sq,er,eq,params,&kernel_score);
  context.phase="stitch_gap_traceback";
  StitchGapResult out;out.engine=dp?"popoa_dp":"ggwfa";out.dp_cells=dp?cells:0;
  OutputLease lease{context};
  std::uint64_t ri=0,qi=0;
  for(const auto& pair:alignment) {
    cg_budget::tick();
    const bool rc=pair.node_id1!=cv::AlignedPair::gap,qc=pair.node_id2!=cv::AlignedPair::gap;
    if((!rc&&!qc)||(rc&&(pair.node_id1!=ri||ri>=r.size()))||(qc&&(pair.node_id2!=qi||qi>=q.size())))
      throw AlignmentError("centromere stitch kernel returned invalid traceback coordinates");
    const char op=!rc?'I':!qc?'D':EqualBase(r[ri],q[qi])?'=':'X';
    Append(out.cigar,op,lease);if(rc)++ri;if(qc)++qi;
  }
  if(ri!=r.size()||qi!=q.size())throw AlignmentError("centromere stitch kernel traceback omitted input bases");
  context.check();
  context.phase="stitch_gap_rescore";
  out.score=ScoreCentromereStitchCigar(r,q,out.cigar,options.score_model);
  if(out.score!=kernel_score)throw AlignmentError("centromere stitch kernel score differs from independent CIGAR rescore");
  context.check();out.peak_kernel_bytes=context.peak;out.wfa_states=context.states;
  out.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-context.start).count();
  return out;
}
} // namespace
#endif

StitchGapResult AlignCentromereStitchGap(std::string_view r,std::string_view q,
                                       const CentromereStitchOptions& options) {
#if !RAMAG_HAVE_CENTROMERE
  (void)r;(void)q;(void)options;
  throw DependencyError("centromere stitching backend was disabled at build time");
#else
  GetScores(options.score_model);
  if(!options.max_dp_cells||!options.gap_memory_limit_bytes||!std::isfinite(options.gap_seconds_limit)||options.gap_seconds_limit<=0)
    throw AlignmentError("invalid centromere stitch gap resource limits");
  if(r.empty()||q.empty())throw AlignmentError("two-sided stitch kernel requires nonempty inputs");
  if(q.size()>std::numeric_limits<std::uint64_t>::max()-2||r.size()>std::numeric_limits<std::uint64_t>::max()-q.size()-2)
    throw AlignmentError("centromere stitch input length overflow");
  const auto total=static_cast<std::uint64_t>(r.size()+q.size());
  const auto rows=static_cast<std::uint64_t>(r.size()+1),cols=static_cast<std::uint64_t>(q.size()+1);
  const auto cells=rows>std::numeric_limits<std::uint64_t>::max()/cols
    ?std::numeric_limits<std::uint64_t>::max():rows*cols;
  const bool dp=cells<=options.max_dp_cells;
  const auto scoring=GetScores(options.score_model);
  const auto unit_bound=std::max({scoring.match,scoring.mismatch,scoring.extend[0]});
  const auto open_bound=*std::max_element(scoring.open.begin(),scoring.open.end());
  // A prefix can always be represented by at most two runs of gaps. Every
  // reachable I/D state has a path bounded by unit_bound*total+2*open_bound.
  // Only the matrix DP uses int32; GGWFA stores its transformed cost in int64.
  const auto scalar_bound=dp?std::int64_t(std::numeric_limits<std::int32_t>::max()/2)
                            :std::numeric_limits<std::int64_t>::max()/4;
  if(total>static_cast<std::uint64_t>((scalar_bound-2*open_bound)/unit_bound)) {
    StitchGapResourceLimit error("kernel_numerical_size_limit");error.engine=dp?"popoa_dp":"ggwfa";
    error.dp_cells=dp?cells:0;throw error;
  }
  cg_budget::Context context;
  context.limit=options.gap_memory_limit_bytes;context.seconds=options.gap_seconds_limit;
  context.callback=[](void* data){
    const auto& opts=*static_cast<const CentromereStitchOptions*>(data);
    if(opts.interruption_callback)opts.interruption_callback(cg_budget::active?cg_budget::active->phase:"stitch_gap_kernel");
  };
  context.callback_data=const_cast<CentromereStitchOptions*>(&options);
  cg_budget::Scope scope(context);
  try {
    context.check();
    auto out=options.score_model==CentromereStitchScoreModel::Centrolign3Piece
      ?Run<3>(r,q,options,context,dp,cells):Run<2>(r,q,options,context,dp,cells);
    assert(context.live==0);
    return out;
  } catch(const cg_budget::Limit& e) {
    assert(context.live==0);
    StitchGapResourceLimit error(e.reason);error.engine=dp?"popoa_dp":"ggwfa";
    error.dp_cells=dp?cells:0;error.peak_kernel_bytes=context.peak;error.wfa_states=context.states;
    error.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-context.start).count();
    throw error;
  } catch(StitchGapResourceLimit& error) {
    // Stage deadline supplied by the assembler. Preserve the exception and
    // reason, adding the work already completed by this internal gap only.
    error.engine=dp?"popoa_dp":"ggwfa";error.dp_cells=dp?cells:0;
    error.peak_kernel_bytes=context.peak;error.wfa_states=context.states;
    error.elapsed_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-context.start).count();
    throw;
  }
#endif
}
} // namespace ramag::detail

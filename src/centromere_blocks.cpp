#include "ramag/centromere_blocks.hpp"
#include "centromere_detail.hpp"
#include "centromere_orientation_region_review.hpp"
#include "centromere_orientation_unified_review.hpp"
#include "centromere_orientation_known_counterchain.hpp"
#include "centromere_orientation_local_direction_repair.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string_view>
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag {
bool CentromereCliBlocks() noexcept {return RAMAG_CENTROMERE_BLOCK_ALIGNMENT!=0;}
CentromereBlockOptions CentromereCliBlockOptions() {
  CentromereBlockOptions o;
#if RAMAG_CENTROMERE_JOINT_BLOCKS
  o.planning_strategy=CentromereBlockStrategy::JointLocalBlocksV2;
#endif
#if RAMAG_CENTROMERE_LOCAL_FRONTIER
  o.planning_strategy=CentromereBlockStrategy::LocalCoherentFrontierV3;
#endif
#if RAMAG_CENTROMERE_ORIENTATION_EVIDENCE
  o.planning_strategy=CentromereBlockStrategy::OrientationEvidenceV9;
#elif RAMAG_CENTROMERE_COUNT_CLUSTER
  o.planning_strategy=CentromereBlockStrategy::CountClusterDirectionsV8;
#elif RAMAG_CENTROMERE_FAST_DIRECTION
  o.planning_strategy=CentromereBlockStrategy::FastLocalDirectionV6;
#elif RAMAG_CENTROMERE_REGIONS_ONLY
  o.planning_strategy=CentromereBlockStrategy::DirectionRegionsOnlyV5;
#elif RAMAG_CENTROMERE_HIERARCHICAL
  o.planning_strategy=CentromereBlockStrategy::HierarchicalDirectionRegionsV4;
#endif
  o.planner_memory_limit_bytes=RAMAG_CENTROMERE_BLOCK_PLANNER_BYTES;
  o.wfa_memory_limit_bytes=RAMAG_CENTROMERE_BLOCK_WFA_BYTES;
  o.wfa_seconds_limit=RAMAG_CENTROMERE_BLOCK_WFA_SECONDS;
  o.alignment_seconds_limit=RAMAG_CENTROMERE_BLOCK_STAGE_SECONDS;
  o.count_window=RAMAG_CENTROMERE_COUNT_WINDOW;o.count_step=RAMAG_CENTROMERE_COUNT_STEP;
  o.count_min_anchors=RAMAG_CENTROMERE_COUNT_MIN_ANCHORS;o.count_dominance_ratio=RAMAG_CENTROMERE_COUNT_DOMINANCE_RATIO;
  o.count_cluster_gap=RAMAG_CENTROMERE_COUNT_CLUSTER_GAP;o.count_inferred_gap=RAMAG_CENTROMERE_COUNT_INFERRED_GAP;
  o.orientation_connected_gaps=RAMAG_CENTROMERE_ORIENTATION_CONNECTED_GAPS!=0;
  o.orientation_bridge_unknown_gaps=RAMAG_CENTROMERE_ORIENTATION_BRIDGE_UNKNOWN_GAPS!=0;
  o.orientation_audit_known_blocks=RAMAG_CENTROMERE_ORIENTATION_AUDIT_KNOWN_BLOCKS!=0;
  o.recover_long_gap_ends=RAMAG_CENTROMERE_RECOVER_LONG_GAP_ENDS!=0;
  o.recover_sparse_pivots=RAMAG_CENTROMERE_RECOVER_SPARSE_PIVOTS!=0;
  o.recover_short_gap_ends=RAMAG_CENTROMERE_RECOVER_SHORT_GAP_ENDS!=0;
  o.recover_mixed_gap_ends=RAMAG_CENTROMERE_RECOVER_MIXED_GAP_ENDS!=0;
  o.recover_sparse_gap_ends=RAMAG_CENTROMERE_RECOVER_SPARSE_GAP_ENDS!=0;
  o.prefer_stronger_sparse_gap_end=RAMAG_CENTROMERE_PREFER_STRONGER_SPARSE_GAP_END!=0;
  o.unify_real_anchor_ends=RAMAG_CENTROMERE_UNIFY_REAL_ANCHOR_ENDS!=0;
  o.review_local_competitive_gaps=RAMAG_CENTROMERE_REVIEW_LOCAL_COMPETITIVE_GAPS!=0;
  o.review_isolated_opposite_anchors=RAMAG_CENTROMERE_REVIEW_ISOLATED_OPPOSITE_ANCHORS!=0;
  o.review_collinear_fragment_connections=RAMAG_CENTROMERE_REVIEW_COLLINEAR_FRAGMENT_CONNECTIONS!=0;
  o.review_collinear_window_gaps=RAMAG_CENTROMERE_REVIEW_COLLINEAR_WINDOW_GAPS!=0;
  o.review_alternate_original_endpoints=RAMAG_CENTROMERE_REVIEW_ALTERNATE_ORIGINAL_ENDPOINTS!=0;
  o.review_neutral_window_gaps=RAMAG_CENTROMERE_REVIEW_NEUTRAL_WINDOW_GAPS!=0;
  o.review_region_chains=RAMAG_CENTROMERE_REVIEW_REGION_CHAINS!=0;
  o.supplement_region_anchors=RAMAG_CENTROMERE_SUPPLEMENT_REGION_ANCHORS!=0;
  o.review_unified_gaps=RAMAG_CENTROMERE_REVIEW_UNIFIED_GAPS!=0;
  o.review_gap_interior_opposition=RAMAG_CENTROMERE_REVIEW_GAP_INTERIOR_OPPOSITION!=0;
  o.review_known_block_counterchains=RAMAG_CENTROMERE_REVIEW_KNOWN_BLOCK_COUNTERCHAINS!=0;
  o.repair_inherited_directions=RAMAG_CENTROMERE_REPAIR_INHERITED_DIRECTIONS!=0;
  o.repair_drift_consistent_gaps=RAMAG_CENTROMERE_REPAIR_DRIFT_CONSISTENT_GAPS!=0;
  o.unified_search_query_limit=RAMAG_CENTROMERE_UNIFIED_QUERY_LIMIT;
  o.unified_search_reference_limit=RAMAG_CENTROMERE_UNIFIED_REFERENCE_LIMIT;
  o.unified_search_input_budget=RAMAG_CENTROMERE_UNIFIED_INPUT_BUDGET;
  return o;
}
std::string_view CentromereBlockStrategyName(CentromereBlockStrategy strategy) noexcept {
  switch(strategy) {
    case CentromereBlockStrategy::InitialDirectionBlocksV1:return "initial-direction-blocks-v1";
    case CentromereBlockStrategy::JointLocalBlocksV2:return "joint-local-blocks-v2";
    case CentromereBlockStrategy::LocalCoherentFrontierV3:return "local-coherent-frontier-v3";
    case CentromereBlockStrategy::HierarchicalDirectionRegionsV4:return "hierarchical-direction-regions-v4";
    case CentromereBlockStrategy::DirectionRegionsOnlyV5:return "direction-regions-only-v5";
    case CentromereBlockStrategy::FastLocalDirectionV6:return "fast-local-direction-v6";
    case CentromereBlockStrategy::CountClusterDirectionsV8:return "count-cluster-directions-v8";
    case CentromereBlockStrategy::OrientationEvidenceV9:return "orientation-evidence-v9";
  }
  return "unknown";
}
namespace {
using Clock=std::chrono::steady_clock;
double Elapsed(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
bool Canonical(char c) {return c=='A'||c=='C'||c=='G'||c=='T';}
char Complement(char c) {
  switch(c) {case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';case 'N':return 'N';}
  throw CentromereInputError("centromere input must be normalized A/C/G/T/N");
}
std::string Reverse(std::string_view text) {
  std::string result(text.size(),'N');
  for(std::size_t i=0;i<text.size();++i) result[text.size()-i-1]=Complement(text[i]);
  return result;
}
void Check(const CentromereAlignmentOptions& o,std::string_view stage) {
  if(o.interruption_callback)o.interruption_callback(stage);
}
void Inputs(const SequenceRecord& r,const SequenceRecord& q,const CentromereAlignmentOptions& o) {
  if(!o.worker_threads || r.bases.empty() || q.bases.empty())
    throw CentromereInputError("centromere sequences and thread budget must be nonempty");
  constexpr Length limit=(INT_MAX-64)/3;
  if(r.size()>limit || q.size()>limit-r.size())throw CentromereInputError("centromere pair exceeds pinned coordinate range");
  for(const auto* sequence:{&r,&q})for(char c:sequence->bases)(void)Complement(c);
}
Position U(const CentromereAnchor& a,Strand strand,Length length) {
  return strand==Strand::Forward ? a.query_begin : length-a.query_begin-a.length;
}
void Chain(std::span<const CentromereAnchor> chain,Strand strand,Length rlength,Length qlength) {
  Position re=0,qe=0;
  for(const auto& a:chain) {
    if(!a.length || a.reference_begin>rlength || a.length>rlength-a.reference_begin ||
       a.query_begin>qlength || a.length>qlength-a.query_begin || !std::isfinite(a.weight) || a.weight<0)
      throw AlignmentError("invalid initial centromere anchor");
    const auto u=U(a,strand,qlength);
    if(a.reference_begin<re || u<qe)throw AlignmentError("initial centromere chain overlaps or is unordered");
    re=a.reference_begin+a.length;qe=u+a.length;
  }
}
Length Support(std::span<const CentromereAnchor> anchors,Strand strand,
               const CentromereBlock& block,Length qlength) {
  const auto qb=strand==Strand::Forward ? block.query_begin : qlength-block.query_end;
  const auto qe=strand==Strand::Forward ? block.query_end : qlength-block.query_begin;
  Length sum=0;
  for(const auto& a:anchors) {
    const auto u=U(a,strand,qlength);
    if(a.reference_begin>=block.reference_end || u>=qe)continue;
    const auto low=std::max(block.reference_begin>a.reference_begin ? block.reference_begin-a.reference_begin : 0,
                            qb>u ? qb-u : 0);
    const auto high=std::min({a.length,block.reference_end-a.reference_begin,qe-u});
    if(low<high)sum+=high-low;
  }
  return sum; // Each source chain is already collinear and nonoverlapping.
}
void Append(std::vector<CigarOp>& cigar,char operation,Length n) {
  if(!n)return;
  if(!cigar.empty() && cigar.back().operation==operation)cigar.back().length+=n;
  else cigar.push_back({operation,n});
}
bool Paired(const AlignmentRecord& a) {
  return std::any_of(a.cigar.begin(),a.cigar.end(),[](const auto& c){return c.operation=='='||c.operation=='X';});
}
void Validate(AlignmentRecord& a,const SequenceRecord& r,const SequenceRecord& q) {
  if(a.reference_begin>a.reference_end || a.reference_end>r.size() ||
     a.query_begin>a.query_end || a.query_end>q.size())throw AlignmentError("block record outside original sequence");
  Length ri=0,qi=0,nm=0;
  std::int64_t cost=0;
  for(const auto& c:a.cigar) {
    if(!c.length)throw AlignmentError("zero block CIGAR run");
    if(c.operation=='=' || c.operation=='X') {
      for(Length i=0;i<c.length;++i) {
        if(ri>=a.reference_end-a.reference_begin || qi>=a.query_end-a.query_begin)
          throw AlignmentError("block CIGAR paired column outside assigned interval");
        const char rb=r.bases[a.reference_begin+ri];
        const char qb=a.strand==Strand::Forward ? q.bases[a.query_begin+qi] : Complement(q.bases[a.query_end-qi-1]);
        if((Canonical(rb)&&rb==qb)!=(c.operation=='='))throw AlignmentError("block CIGAR disagrees with original bases at reference="+
          std::to_string(a.reference_begin+ri)+" query_offset="+std::to_string(qi)+" strand="+
          (a.strand==Strand::Forward ? "+" : "-")+" reference_base="+rb+" oriented_query_base="+qb);
        ++ri;++qi;
      }
      if(c.operation=='X'){cost+=3*static_cast<std::int64_t>(c.length);nm+=c.length;}
    } else if(c.operation=='I' || c.operation=='D') {
      const auto n=static_cast<std::int64_t>(c.length);
      cost+=std::min(4+2*n,12+n);nm+=c.length;
      if(c.operation=='I')qi+=c.length;else ri+=c.length;
    } else throw AlignmentError("invalid block CIGAR operation");
    if(ri>a.reference_end-a.reference_begin || qi>a.query_end-a.query_begin)
      throw AlignmentError("block CIGAR overconsumes its interval");
  }
  if(ri!=a.reference_end-a.reference_begin || qi!=a.query_end-a.query_begin)
    throw AlignmentError("block CIGAR consumption does not close");
  a.score=-cost;a.edit_distance=nm;
}
AlignmentRecord Piece(const CentromereBlock& b,Position rb,Position re,Position qb,Position qe,
                      std::vector<CigarOp> cigar,const SequenceRecord& r,const SequenceRecord& q) {
  AlignmentRecord a;
  a.reference_id=r.numeric_id;a.query_id=q.numeric_id;a.strand=b.strand;
  a.reference_begin=b.reference_begin+rb;a.reference_end=b.reference_begin+re;
  a.query_begin=b.strand==Strand::Forward ? b.query_begin+qb : b.query_end-qe;
  a.query_end=b.strand==Strand::Forward ? b.query_begin+qe : b.query_end-qb;
  a.cigar=std::move(cigar);Validate(a,r,q);return a;
}
bool Adjacent(const AlignmentRecord& a,const AlignmentRecord& b) {
  return a.strand==b.strand && a.reference_end==b.reference_begin &&
    (a.strand==Strand::Forward ? a.query_end==b.query_begin : a.query_begin==b.query_end);
}
void Missing(CentromereBlockResult& out,const CentromereBlock& parent,
             Position rb,Position re,Position qb,Position qe,std::string reason,bool unpaired=false) {
  auto b=parent;b.reference_begin=parent.reference_begin+rb;b.reference_end=parent.reference_begin+re;
  b.query_begin=parent.strand==Strand::Forward ? parent.query_begin+qb : parent.query_end-qe;
  b.query_end=parent.strand==Strand::Forward ? parent.query_begin+qe : parent.query_end-qb;
  b.status=unpaired ? "unpaired" : "unresolved";b.reason=std::move(reason);
  out.unresolved.push_back({std::move(b),unpaired});
  if(!unpaired)out.complete=false;
}
void Ledger(const CentromereBlockResult& out,Length reference_length,Length query_length) {
  for(bool reference:{true,false}) {
    std::vector<std::pair<Position,Position>> intervals;
    for(const auto& a:out.alignment.alignments) intervals.emplace_back(reference ? a.reference_begin : a.query_begin,
                                                                     reference ? a.reference_end : a.query_end);
    for(const auto& missing:out.unresolved) {
      const auto& a=missing.interval;
      const auto b=reference ? a.reference_begin : a.query_begin,e=reference ? a.reference_end : a.query_end;
      if(b<e)intervals.emplace_back(b,e);
    }
    std::sort(intervals.begin(),intervals.end());Position end=0;
    for(const auto& interval:intervals) {
      if(interval.first!=end || interval.second<interval.first)throw AlignmentError("block ledger has a hole or overlapping ownership");
      end=interval.second;
    }
    if(end!=(reference ? reference_length : query_length))throw AlignmentError("block ledger does not explain entire input");
  }
}
} // namespace

namespace {
std::vector<CentromereAnchor> InitialEvidence(const SequenceRecord& r,
    const SequenceRecord& q,const CentromereAlignmentOptions& o,
    std::vector<CentromerePlanningAnchor>* candidates=nullptr,Strand strand=Strand::Forward,
    std::uint64_t memory_limit=256000000) {
  Inputs(r,q,o);
#if RAMAG_HAVE_CENTROMERE
  Check(o,"centromere-initial");
  std::vector<rama_source::SequenceInfo> sequences{{r.bases,r.name,static_cast<uint_t>(r.size())},
                                                 {q.bases,q.name,static_cast<uint_t>(q.size())}};
  rama_source::AnchorFinder finder(sequences,"",o.worker_threads==1 ? 0 : o.worker_threads,false,false,100,
                                   [&]{Check(o,"centromere-initial");});
  rama_source::RareCandidateSink observer;
  if(candidates)observer=[&](const rama_source::RareCandidate& c) {
    const auto& a=c.pair;
    if(!c.first_hit_level || a.second_pos<r.size()+1)throw AlignmentError("invalid root evidence source");
    const auto u=a.second_pos-r.size()-1;
    if(a.first_pos>r.size() || a.match_length>r.size()-a.first_pos || u>q.size() || a.match_length>q.size()-u)
      throw AlignmentError("root candidate outside input");
    for(Length i=0;i<a.match_length;++i)
      if(!Canonical(r.bases[a.first_pos+i]) || r.bases[a.first_pos+i]!=q.bases[u+i])
        throw AlignmentError("root candidate is not canonical exact");
    if(candidates->size()==candidates->capacity()) {
      const auto capacity=std::max<std::size_t>(1024,candidates->capacity()*2);
      if(capacity>memory_limit/sizeof(CentromerePlanningAnchor) || candidates->capacity()>
         memory_limit/sizeof(CentromerePlanningAnchor)-capacity)
        throw AlignmentError("root candidate capacity budget exceeded; no truncation");
      candidates->reserve(capacity);
    }
    candidates->push_back({{a.first_pos,strand==Strand::Forward ? u : q.size()-u-a.match_length,
      a.match_length,a.weight},strand,static_cast<std::uint32_t>(c.reference_count),
      static_cast<std::uint32_t>(c.query_count),static_cast<std::uint32_t>(c.level)});
  };
  const auto pairs=finder.findInitialAnchors(observer);
  std::vector<CentromereAnchor> result;result.reserve(pairs.size());
  for(const auto& a:pairs) {
    if(a.second_pos<r.size()+1)throw AlignmentError("invalid initial query offset");
    const auto qb=a.second_pos-r.size()-1;
    if(a.first_pos>r.size() || a.match_length>r.size()-a.first_pos || qb>q.size() || a.match_length>q.size()-qb)
      throw AlignmentError("initial match outside input");
    for(Length i=0;i<a.match_length;++i)
      if(!Canonical(r.bases[a.first_pos+i]) || r.bases[a.first_pos+i]!=q.bases[qb+i])
        throw AlignmentError("initial block anchor is not exact");
    result.push_back({a.first_pos,qb,a.match_length,a.weight});
  }
  Chain(result,Strand::Forward,r.size(),q.size());return result;
#else
  (void)candidates;(void)strand;(void)memory_limit;
  throw DependencyError("centromere block backend requires Linux source RaMA backend");
#endif
}
} // namespace
std::vector<CentromereAnchor> FindInitialCentromereAnchors(const SequenceRecord& r,
    const SequenceRecord& q,const CentromereAlignmentOptions& o) {
  return InitialEvidence(r,q,o);
}
std::vector<CentromerePlanningAnchor> detail::SearchCentromereRegionCandidates(
    const SequenceRecord& reference,const SequenceRecord& query,const CentromereRegionSearch& task,
    Strand strand,const CentromereBlockOptions& options) {
  if(task.reference_begin>=task.reference_end || task.reference_end>reference.size() ||
     task.query_begin>=task.query_end || task.query_end>query.size())
    throw CentromereInputError("local source-RaMA rectangle outside input");
  SequenceRecord r{reference.numeric_id,reference.name,reference.header,
    reference.bases.substr(task.reference_begin,task.reference_end-task.reference_begin)};
  SequenceRecord q{query.numeric_id,query.name,query.header,
    query.bases.substr(task.query_begin,task.query_end-task.query_begin)};
  if(strand==Strand::Reverse)q.bases=Reverse(q.bases);
  std::vector<CentromerePlanningAnchor> local;
  (void)InitialEvidence(r,q,options,&local,strand,options.planner_memory_limit_bytes);
  for(auto& x:local){x.anchor.reference_begin+=task.reference_begin;x.anchor.query_begin+=task.query_begin;}
  return local;
}

CentromereBlockResult PlanCentromereBlocks(Length rlength,Length qlength,
    std::vector<CentromereAnchor> forward,std::vector<CentromereAnchor> reverse) {
  if(!rlength || !qlength)throw CentromereInputError("block planner requires nonempty input");
  Chain(forward,Strand::Forward,rlength,qlength);Chain(reverse,Strand::Reverse,rlength,qlength);
  CentromereBlockResult out;Length fs=0,rs=0;
  for(const auto& a:forward)fs+=a.length;
  for(const auto& a:reverse)rs+=a.length;
  out.background=rs>fs ? Strand::Reverse : Strand::Forward;out.background_tied=rs==fs;
  out.initial_forward=std::move(forward);out.initial_reverse=std::move(reverse);
  const auto& backbone=out.background==Strand::Forward ? out.initial_forward : out.initial_reverse;
  Position rb=0,ub=0;
  for(std::size_t i=0;i<=backbone.size();++i) {
    const auto re=i<backbone.size() ? backbone[i].reference_begin : rlength;
    const auto ue=i<backbone.size() ? U(backbone[i],out.background,qlength) : qlength;
    if(re<rb || ue<ub)throw AlignmentError("background source chain is not monotonic");
    if(rb<re || ub<ue) {
      CentromereBlock b;b.id=out.blocks.size()+1;b.reference_begin=rb;b.reference_end=re;
      b.query_begin=out.background==Strand::Forward ? ub : qlength-ue;
      b.query_end=out.background==Strand::Forward ? ue : qlength-ub;
      b.forward_support=Support(out.initial_forward,Strand::Forward,b,qlength);
      b.reverse_support=Support(out.initial_reverse,Strand::Reverse,b,qlength);
      b.tentative=b.forward_support==b.reverse_support;
      b.strand=b.tentative ? out.background : b.reverse_support>b.forward_support ? Strand::Reverse : Strand::Forward;
      b.inversion_candidate=i>0 && i<backbone.size() && b.strand!=out.background;
      b.reason=b.tentative ? "background_fallback" : b.inversion_candidate ? "opposite_support_between_background_anchors" : "initial_direction_support";
      out.blocks.push_back(std::move(b));
    }
    if(i==backbone.size())break;
    const auto& a=backbone[i];
    CentromereBlock b;b.id=out.blocks.size()+1;b.skeleton=true;b.strand=out.background;
    b.reference_begin=a.reference_begin;b.reference_end=a.reference_begin+a.length;
    b.query_begin=a.query_begin;b.query_end=a.query_begin+a.length;
    b.forward_support=out.background==Strand::Forward ? a.length : 0;
    b.reverse_support=out.background==Strand::Reverse ? a.length : 0;b.reason="exact_background_anchor";
    out.blocks.push_back(std::move(b));rb=a.reference_begin+a.length;ub=U(a,out.background,qlength)+a.length;
  }
  return out;
}

CentromereBlockResult FindCentromereDirectionRegions(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereBlockOptions& options) {
  if(options.repair_drift_consistent_gaps && !options.repair_inherited_directions)
    throw CentromereInputError("drift-consistent repair requires inherited-direction repair");
  if(options.repair_inherited_directions && !options.review_known_block_counterchains)
    throw CentromereInputError("inherited-direction repair requires the frozen known-counterchain baseline");
  if(options.review_gap_interior_opposition && !options.review_unified_gaps)
    throw CentromereInputError("gap-interior opposition review requires unified gap review");
  if(options.review_unified_gaps && (!options.review_region_chains || !options.supplement_region_anchors))
    throw CentromereInputError("unified gap review requires the frozen complete-region baseline");
  if(options.supplement_region_anchors && !options.review_region_chains)
    throw CentromereInputError("local region search requires complete region review");
  if(options.review_region_chains && !options.review_neutral_window_gaps)
    throw CentromereInputError("region review requires the frozen neutral-window baseline");
  if(options.prefer_stronger_sparse_gap_end && !options.recover_sparse_gap_ends)
    throw CentromereInputError("stronger sparse-gap endpoint preference requires sparse-gap endpoint recovery");
  if(options.review_isolated_opposite_anchors && !options.review_local_competitive_gaps)
    throw CentromereInputError("isolated opposite review requires frozen local review layers");
  if(options.review_collinear_window_gaps && !options.review_collinear_fragment_connections)
    throw AlignmentError("collinear window-gap review requires collinear fragment review");
  if(options.review_alternate_original_endpoints && !options.review_collinear_window_gaps)
    throw CentromereInputError("alternate endpoint review requires frozen collinear window-gap review");
  if(options.review_neutral_window_gaps && !options.review_alternate_original_endpoints)
    throw CentromereInputError("neutral window review requires bounded alternate original endpoint review");
  if(options.review_collinear_fragment_connections && !options.review_isolated_opposite_anchors)
    throw CentromereInputError("collinear fragment review requires frozen isolated review layers");
  if(options.review_local_competitive_gaps && !options.unify_real_anchor_ends)
    throw CentromereInputError("local competitive review requires frozen unified endpoints");
  if(options.unify_real_anchor_ends && (!options.prefer_stronger_sparse_gap_end ||
     !options.orientation_bridge_unknown_gaps || !options.orientation_audit_known_blocks ||
     options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9))
    throw CentromereInputError("unified real endpoints require the frozen v9q layers");
  if(options.recover_short_gap_ends && !options.recover_sparse_pivots)
    throw CentromereInputError("short-gap endpoint recovery requires the sparse-pivot orientation stack");
  if(options.recover_mixed_gap_ends && !options.recover_short_gap_ends)
    throw CentromereInputError("mixed-gap endpoint recovery requires the short-gap orientation stack");
  if(options.recover_sparse_gap_ends && !options.recover_mixed_gap_ends)
    throw CentromereInputError("sparse-gap endpoint recovery requires the mixed-gap orientation stack");
  if(options.recover_sparse_pivots &&
     (options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9 ||
      !options.orientation_connected_gaps || !options.orientation_bridge_unknown_gaps ||
      !options.orientation_audit_known_blocks || !options.recover_long_gap_ends))
    throw CentromereInputError("sparse root pivots require connected, unknown-gap, audited and long-gap orientation-evidence-v9");
  if(options.orientation_audit_known_blocks && options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9)
    throw CentromereInputError("known-block orientation audit requires orientation-evidence-v9");
  if(options.recover_long_gap_ends &&
     (!options.orientation_connected_gaps || options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9))
    throw CentromereInputError("long-gap endpoint recovery requires connected orientation-evidence-v9");
  if(options.orientation_bridge_unknown_gaps &&
     (!options.orientation_connected_gaps || options.planning_strategy!=CentromereBlockStrategy::OrientationEvidenceV9))
    throw CentromereInputError("unknown-gap bridging requires connected orientation-evidence-v9");
  Inputs(reference,query,options);
  const auto start=Clock::now();
  const bool joint=options.planning_strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1;
  std::vector<CentromerePlanningAnchor> candidates;
  if(options.progress_callback)options.progress_callback("centromere-initial-forward",0,2);
  auto forward=InitialEvidence(reference,query,options,joint ? &candidates : nullptr,
    Strand::Forward,options.planner_memory_limit_bytes);
  std::vector<CentromereAnchor> reverse;
  {
    SequenceRecord reversed{query.numeric_id,query.name,query.header,Reverse(query.bases)};
    if(options.progress_callback)options.progress_callback("centromere-initial-reverse",1,2);
    reverse=InitialEvidence(reference,reversed,options,joint ? &candidates : nullptr,
      Strand::Reverse,options.planner_memory_limit_bytes);
  }
  for(auto& a:reverse)a.query_begin=query.size()-a.query_begin-a.length;
  const double initial_seconds=Elapsed(start);const auto planning=Clock::now();
  CentromereBlockResult out;
  if(joint) {
    if(options.progress_callback)options.progress_callback("centromere-joint-planning",0,1);
    Length fs=0,rs=0;for(const auto& a:forward)fs+=a.length;for(const auto& a:reverse)rs+=a.length;
    CentromereBlockResult discovered;discovered.strategy=options.planning_strategy;
    discovered.root_candidates=std::move(candidates);discovered.initial_forward=std::move(forward);
    discovered.initial_reverse=std::move(reverse);
    if(options.observation_callback)options.observation_callback(discovered);
    if(options.planning_strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4 ||
       options.planning_strategy==CentromereBlockStrategy::DirectionRegionsOnlyV5) {
      try {
        out=detail::ReplanDirectionRegions(reference.size(),query.size(),discovered.root_candidates,options,
          rs>fs ? Strand::Reverse : Strand::Forward);
      } catch(const detail::CentromereBudgetExceeded& e) {
        out.refinement_partial_reason=e.what();
      } catch(const std::bad_alloc&) {out.refinement_partial_reason="initial_planner_allocation_failed";}
      if(out.refinement_partial_reason.empty() && options.planning_strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4)
        out=detail::RefineDirectionRegions(reference,query,std::move(out),options,rs>fs ? Strand::Reverse : Strand::Forward);
      else if(!out.refinement_partial_reason.empty()) {
        out.strategy=options.planning_strategy;out.root_candidates=std::move(discovered.root_candidates);
        out.root_candidate_count=out.root_candidates.size();
        CentromereBlock unknown;unknown.id=1;unknown.reference_end=reference.size();unknown.query_end=query.size();
        unknown.direction_known=false;unknown.tentative=true;unknown.status="unresolved";unknown.reason="direction_unresolved";
        out.blocks.push_back(std::move(unknown));
      }
    } else out=detail::ReplanDirectionRegions(reference.size(),query.size(),discovered.root_candidates,options,
      rs>fs ? Strand::Reverse : Strand::Forward,reference.bases,query.bases);
    out.initial_forward=std::move(discovered.initial_forward);out.initial_reverse=std::move(discovered.initial_reverse);
  } else out=PlanCentromereBlocks(reference.size(),query.size(),std::move(forward),std::move(reverse));
  auto local_search=[&](const CentromereRegionSearch& task,Strand strand) {
    return detail::SearchCentromereRegionCandidates(reference,query,task,strand,options);
  };
  if(options.review_region_chains)detail::ReviewOrientationRegions(reference.bases,query.bases,out,options,local_search);
  if(options.review_unified_gaps)detail::ReviewUnifiedOrientationGaps(reference.bases,query.bases,out,options,local_search);
  if(options.review_known_block_counterchains)detail::ReviewKnownBlockCounterchains(query.bases,out,options);
  if(options.repair_inherited_directions)detail::RepairInheritedDirections(query.bases,out,options);
  out.initial_seconds=initial_seconds;out.planning_seconds=Elapsed(planning);
  if(options.observation_callback)options.observation_callback(out);
  return out;
}
CentromereBlockResult AlignCentromeresWithDirections(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereBlockOptions& options) {
  if(options.planning_strategy==CentromereBlockStrategy::CountClusterDirectionsV8 ||
     options.planning_strategy==CentromereBlockStrategy::OrientationEvidenceV9)
    throw CentromereInputError("count-cluster-directions-v8 is direction-only; use FindCentromereDirectionRegions");
  if(!options.wfa_memory_limit_bytes || !std::isfinite(options.wfa_seconds_limit) || options.wfa_seconds_limit<0 ||
     !std::isfinite(options.alignment_seconds_limit) || options.alignment_seconds_limit<0)
    throw CentromereInputError("invalid centromere block budgets");
  const auto start=Clock::now();
  auto out=FindCentromereDirectionRegions(reference,query,options);
  const bool joint=options.planning_strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1;
  const auto alignment_start=Clock::now();
  auto stage_check=[&](std::string_view phase) {
    Check(options,phase);
    if(Elapsed(alignment_start)>=options.alignment_seconds_limit)
      throw detail::CentromereBudgetExceeded("alignment_stage_time_budget");
  };
  std::vector<AlignmentRecord> pieces;
  for(auto& block:out.blocks) {
    const auto piece_begin=pieces.size();
    Check(options,"centromere-blocks");
    if(options.progress_callback)options.progress_callback("centromere-blocks",block.id-1,out.blocks.size());
    const auto rn=block.reference_end-block.reference_begin,qn=block.query_end-block.query_begin;
    if(!block.direction_known) {
      out.unresolved.push_back({block,false});out.complete=false;continue;
    }
    if(block.skeleton) {
      pieces.push_back(Piece(block,0,rn,0,qn,{{'=',rn}},reference,query));
      out.refined.push_back({{block.reference_begin,block.query_begin,rn,static_cast<double>(rn)},block.strand,block.id});
      if(block.region_id)out.regions.at(block.region_id-1).paired_columns+=rn;
      block.status="complete";continue;
    }
    std::vector<CentromereAnchor> anchors;
    struct FailedGap {Position rb,re,qb,qe;std::string reason;};
    std::vector<FailedGap> failed_recursions;
    SequenceRecord r{reference.numeric_id,reference.name,reference.header,reference.bases.substr(block.reference_begin,rn)};
    SequenceRecord q{query.numeric_id,query.name,query.header,query.bases.substr(block.query_begin,qn)};
    if(block.strand==Strand::Reverse)q.bases=Reverse(q.bases);
    const auto recursion_start=Clock::now();
    try {
      if(!joint)stage_check("centromere-block-recursion");
      if(rn && qn) {
        CentromereAlignmentOptions local=options;local.progress_callback={};
        local.interruption_callback=stage_check;
        if(!joint)anchors=FindCentromereAnchors(r,q,local).final;
        else {
          // The selected planning chain is immutable. Recurse only in its
          // gaps, never reselect a different correspondence across its anchors.
          Position rb=0,qb=0;
          for(std::size_t i=0;i<=block.planning_anchors.size();++i) {
            const auto re=i<block.planning_anchors.size() ? block.planning_anchors[i].reference_begin-block.reference_begin : rn;
            const auto qe=i<block.planning_anchors.size() ?
              (block.strand==Strand::Forward ? block.planning_anchors[i].query_begin-block.query_begin :
                block.query_end-block.planning_anchors[i].query_begin-block.planning_anchors[i].length) : qn;
            if(re<rb || qe<qb)throw AlignmentError("joint planning anchors are not collinear");
            if(re>rb && qe>qb) {
              try {
                stage_check("centromere-joint-gap-recursion");
                SequenceRecord gr{r.numeric_id,r.name,r.header,r.bases.substr(rb,re-rb)};
                SequenceRecord gq{q.numeric_id,q.name,q.header,q.bases.substr(qb,qe-qb)};
                auto found=FindCentromereAnchors(gr,gq,local).final;
                for(auto& a:found){a.reference_begin+=rb;a.query_begin+=qb;anchors.push_back(a);}
              } catch(const detail::CentromereBudgetExceeded& e) {
                failed_recursions.push_back({rb,re,qb,qe,e.what()});
              } catch(const std::bad_alloc&) {
                failed_recursions.push_back({rb,re,qb,qe,"recursion_allocation_failed"});
              }
            }
            if(i<block.planning_anchors.size()) {
              const auto& a=block.planning_anchors[i];anchors.push_back({re,qe,a.length,a.weight});
              rb=re+a.length;qb=qe+a.length;
            }
          }
        }
      }
    } catch(const detail::CentromereBudgetExceeded& e) {
      Missing(out,block,0,rn,0,qn,e.what());block.status="unresolved";block.reason=e.what();
      out.recursion_seconds+=Elapsed(recursion_start);continue;
    } catch(const std::bad_alloc&) {
      Missing(out,block,0,rn,0,qn,"recursion_allocation_failed");block.status="unresolved";block.reason="recursion_allocation_failed";
      out.recursion_seconds+=Elapsed(recursion_start);continue;
    }
    out.recursion_seconds+=Elapsed(recursion_start);
    bool incomplete=false;Position rb=0,qb=0;
    for(std::size_t i=0;i<=anchors.size();++i) {
      const auto re=i<anchors.size() ? anchors[i].reference_begin : rn;
      const auto qe=i<anchors.size() ? anchors[i].query_begin : qn;
      if(re<rb || qe<qb)throw AlignmentError("local source recursion is not monotonic");
      if(re>rb || qe>qb) {
        const auto failed=std::find_if(failed_recursions.begin(),failed_recursions.end(),[&](const auto& x){
          return x.rb==rb && x.re==re && x.qb==qb && x.qe==qe;});
        if(failed!=failed_recursions.end()) {
          Missing(out,block,rb,re,qb,qe,failed->reason);incomplete=true;
        } else {
        const auto gap_start=Clock::now();
        try {
          const bool actual_wfa=re>rb && qe>qb && !((re-rb<=5 && qe-qb>100)||(qe-qb<=5 && re-rb>100));
          CentromereAlignmentOptions local=options;local.progress_callback={};
          local.interruption_callback=[&](std::string_view phase) {
            stage_check(phase);
            if(Elapsed(gap_start)>=options.wfa_seconds_limit)throw detail::CentromereBudgetExceeded("wfa_time_budget");
          };
          if(actual_wfa) {stage_check("centromere-block-wfa");++out.wfa_calls;}
          auto cigar=detail::AlignRamaGap(std::string_view(r.bases).substr(rb,re-rb),
                                        std::string_view(q.bases).substr(qb,qe-qb),local,options.wfa_memory_limit_bytes);
          pieces.push_back(Piece(block,rb,re,qb,qe,std::move(cigar),reference,query));
        } catch(const detail::CentromereBudgetExceeded& e) {
          Missing(out,block,rb,re,qb,qe,e.what());incomplete=true;
        } catch(const std::bad_alloc&) {
          Missing(out,block,rb,re,qb,qe,"wfa_allocation_failed");incomplete=true;
        }
        out.wfa_seconds+=Elapsed(gap_start);
        }
      }
      if(i<anchors.size()) {
        const auto& a=anchors[i];
        pieces.push_back(Piece(block,a.reference_begin,a.reference_begin+a.length,a.query_begin,a.query_begin+a.length,
                               {{'=',a.length}},reference,query));
        const auto& piece=pieces.back();
        out.refined.push_back({{piece.reference_begin,piece.query_begin,a.length,a.weight},block.strand,block.id});
        rb=a.reference_begin+a.length;qb=a.query_begin+a.length;
      }
    }
    block.status=incomplete ? "partial" : "complete";
    if(block.region_id) {
      auto& region=out.regions.at(block.region_id-1);
      for(auto i=piece_begin;i<pieces.size();++i)for(const auto& c:pieces[i].cigar) {
        if(c.operation=='I')region.insertions+=c.length;
        else if(c.operation=='D')region.deletions+=c.length;
        else region.paired_columns+=c.length;
      }
    }
    if(options.observation_callback)options.observation_callback(out);
  }
  // Join only truly contiguous consumption. An unresolved gap can never be
  // removed and its two flanks silently stitched together.
  std::vector<AlignmentRecord> joined;
  for(auto& a:pieces) {
    if(!joined.empty() && Adjacent(joined.back(),a)) {
      auto& b=joined.back();b.reference_end=a.reference_end;
      b.query_begin=std::min(b.query_begin,a.query_begin);b.query_end=std::max(b.query_end,a.query_end);
      for(const auto& c:a.cigar)Append(b.cigar,c.operation,c.length);
    } else joined.push_back(std::move(a));
  }
  for(auto& a:joined) {
    Validate(a,reference,query);
    if(a.reference_begin<a.reference_end && a.query_begin<a.query_end && Paired(a))
      out.alignment.alignments.push_back(std::move(a));
    else {
      CentromereBlock b;b.reference_begin=a.reference_begin;b.reference_end=a.reference_end;
      b.query_begin=a.query_begin;b.query_end=a.query_end;b.strand=a.strand;b.status="unpaired";b.reason="no_paired_columns";
      out.unresolved.push_back({std::move(b),true});
    }
  }
  std::size_t primary=0;Length best=0;
  for(std::size_t i=0;i<out.alignment.alignments.size();++i) {
    Length matches=0;
    for(const auto& op:out.alignment.alignments[i].cigar)if(op.operation=='=')matches+=op.length;
    if(matches>best){best=matches;primary=i;}
  }
  if(!out.alignment.alignments.empty())out.alignment.alignments[primary].primary=true;
  Ledger(out,reference.size(),query.size());
  auto& stats=out.alignment.statistics;stats.reference_contigs=stats.query_contigs=1;
  stats.reference_bases=reference.size();stats.query_bases=query.size();stats.alignment_count=out.alignment.alignments.size();
  stats.actual_seed_route="rama-source-initial-both-v1";stats.chaining_route=CentromereBlockStrategyName(out.strategy);
  stats.seed_seconds=out.initial_seconds;stats.chain_and_extension_seconds=out.recursion_seconds+out.wfa_seconds;
  stats.total_seconds=Elapsed(start);stats.seed_requested_threads=options.worker_threads;
  stats.seed_worker_threads=options.worker_threads;stats.extension_worker_threads=1;
  if(options.observation_callback)options.observation_callback(out);
  return out;
}

namespace {
char Sign(Strand strand) {return strand==Strand::Forward ? '+' : '-';}
void Header(std::ostream& o,bool joint) {o<<"block_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tskeleton\ttentative\tinversion_candidate\tforward_support\treverse_support\tstatus\treason";if(joint)o<<"\tdirection_known";o<<'\n';}
void Row(std::ostream& o,const CentromereBlock& b,bool joint) {
  o<<b.id<<'\t'<<b.reference_begin<<'\t'<<b.reference_end<<'\t'<<b.query_begin<<'\t'<<b.query_end<<'\t'<<(b.direction_known ? Sign(b.strand) : '.')<<'\t'
   <<b.skeleton<<'\t'<<b.tentative<<'\t'<<b.inversion_candidate<<'\t'<<b.forward_support<<'\t'<<b.reverse_support<<'\t'<<b.status<<'\t'<<b.reason;if(joint)o<<'\t'<<b.direction_known;o<<'\n';
}
void AnchorHeader(std::ostream& o) {o<<"anchor_id\tblock_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tmatch_length\tweight\n";o<<std::setprecision(std::numeric_limits<double>::max_digits10);}
void AnchorRow(std::ostream& o,const CentromereAnchor& a,Strand strand,std::uint64_t id,std::uint64_t block) {
  o<<id<<'\t'<<block<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\n';
}
}
void WriteCentromereBlocks(std::ostream& o,const CentromereBlockResult& r) {
  const bool joint=r.strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1;
  if(r.strategy!=CentromereBlockStrategy::HierarchicalDirectionRegionsV4 && r.strategy!=CentromereBlockStrategy::DirectionRegionsOnlyV5 &&
     r.strategy!=CentromereBlockStrategy::CountClusterDirectionsV8){Header(o,joint);for(const auto& b:r.blocks)Row(o,b,joint);return;}
  o<<"block_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tskeleton\ttentative\tinversion_candidate\tforward_support\treverse_support\tstatus\treason\tdirection_known\tregion_id\n";
  for(const auto& b:r.blocks) {
    std::ostringstream line;Row(line,b,true);auto text=line.str();text.pop_back();o<<text<<'\t'<<b.region_id<<'\n';
  }
}
void WriteCentromereInitialAnchors(std::ostream& o,const CentromereBlockResult& r) {
  AnchorHeader(o);std::uint64_t id=0;
  for(const auto& a:r.initial_forward)AnchorRow(o,a,Strand::Forward,++id,0);
  for(const auto& a:r.initial_reverse)AnchorRow(o,a,Strand::Reverse,++id,0);
}
void WriteCentromereRefinedAnchors(std::ostream& o,const CentromereBlockResult& r) {
  AnchorHeader(o);std::uint64_t id=0;for(const auto& a:r.refined)AnchorRow(o,a.anchor,a.strand,++id,a.block_id);
}
void WriteCentromereUnresolved(std::ostream& o,const CentromereBlockResult& r) {const bool joint=r.strategy!=CentromereBlockStrategy::InitialDirectionBlocksV1;Header(o,joint);for(const auto& a:r.unresolved)Row(o,a.interval,joint);}
void WriteCentromereRootCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"candidate_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tmatch_length\tweight\treference_count\tquery_count\trare_level\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  const auto size=(r.strategy==CentromereBlockStrategy::HierarchicalDirectionRegionsV4 || r.strategy==CentromereBlockStrategy::DirectionRegionsOnlyV5) ? r.root_candidate_count : r.root_candidates.size();
  for(std::size_t i=0;i<size;++i) {
    const auto& c=r.root_candidates[i];const auto& a=c.anchor;
    o<<i+1<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length
     <<'\t'<<Sign(c.strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\t'<<c.reference_count<<'\t'<<c.query_count<<'\t'<<c.rare_level<<'\n';
  }
}
void WriteCentromereDirectionRegions(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\texact_support\tsecondary_score\tcompeting_mappings\tblock_count\tleft_boundary_width\tright_boundary_width\tpaired_columns\tinsertions\tdeletions\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(const auto& x:r.regions)o<<x.id<<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<Sign(x.strand)<<'\t'
    <<x.support<<'\t'<<x.secondary_score<<'\t'<<x.competing_mappings<<'\t'<<x.block_count<<'\t'<<x.left_boundary_width<<'\t'<<x.right_boundary_width
    <<'\t'<<x.paired_columns<<'\t'<<x.insertions<<'\t'<<x.deletions<<'\n';
}
void WriteCentromereRefinementRegions(std::ostream& o,const CentromereBlockResult& r) {
  o<<"search_region_id\tround\treference_begin\treference_end\tquery_begin\tquery_end\tnew_mappings\tstatus\n";
  for(const auto& x:r.refinement_regions)o<<x.id<<'\t'<<x.round<<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<x.new_mappings<<'\t'<<x.status<<'\n';
}
void WriteCentromereLocalCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"candidate_id\tsearch_region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tmatch_length\tweight\tlocal_reference_count\tlocal_query_count\trare_level\tfull_reference_count\tfull_query_count\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(const auto& x:r.local_candidates) {
    const auto& c=r.root_candidates.at(x.candidate_index);const auto& a=c.anchor;
    o<<x.candidate_index+1<<'\t'<<x.region_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length
      <<'\t'<<Sign(c.strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\t'<<c.reference_count<<'\t'<<c.query_count<<'\t'<<c.rare_level<<'\t'
      <<x.full_reference_count<<'\t'<<x.full_query_count<<'\n';
  }
}
void WriteCentromerePlanningAnchors(std::ostream& o,const CentromereBlockResult& r) {
  AnchorHeader(o);std::uint64_t id=0;
  for(const auto& b:r.blocks)for(const auto& a:b.planning_anchors)AnchorRow(o,a,b.strand,++id,b.id);
}
void WriteCentromereCountWindows(std::ostream& o,const CentromereBlockResult& r) {
  o<<"window_id\tquery_begin\tquery_end\tforward_count\treverse_count\tlabel\n";
  for(const auto& w:r.count_windows)o<<w.id<<'\t'<<w.begin<<'\t'<<w.end<<'\t'<<w.forward<<'\t'<<w.reverse<<'\t'<<w.label<<'\n';
}
void WriteCentromereCountNodes(std::ostream& o,const CentromereBlockResult& r) {
  o<<"node_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tweight\tpredecessor\tfirst\tcount\texact_support\tpath_weight\tneutral\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(std::size_t i=0;i<r.count_nodes.size();++i) {
    const auto& n=r.count_nodes[i];const auto& a=n.anchor;
    o<<i+1<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length
      <<'\t'<<Sign(n.strand)<<'\t'<<a.weight<<'\t'<<n.predecessor<<'\t'<<n.first<<'\t'<<n.count<<'\t'<<n.support<<'\t'<<n.weight<<'\t'<<n.neutral<<'\n';
  }
}
void WriteCentromereCountCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"candidate_id\tfirst\tendpoint\tcount\tforward_count\treverse_count\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\texact_support\tweight\teligible\tselected\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(const auto& c:r.count_candidates)o<<c.id<<'\t'<<c.first<<'\t'<<c.endpoint<<'\t'<<c.count<<'\t'<<c.forward<<'\t'<<c.reverse
    <<'\t'<<c.rb<<'\t'<<c.re<<'\t'<<c.qb<<'\t'<<c.qe<<'\t'<<Sign(c.strand)<<'\t'<<c.support<<'\t'<<c.weight<<'\t'<<c.eligible<<'\t'<<c.selected<<'\n';
}
void WriteCentromereCountGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\tleft_node\tright_node\twindow_id\tforward_count\treverse_count\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\taccepted\treason\n";
  for(const auto& g:r.count_gaps)o<<g.region_id<<'\t'<<g.left_node<<'\t'<<g.right_node<<'\t'<<g.window_id<<'\t'<<g.forward<<'\t'<<g.reverse
    <<'\t'<<g.rb<<'\t'<<g.re<<'\t'<<g.qb<<'\t'<<g.qe<<'\t'<<Sign(g.strand)<<'\t'<<g.accepted<<'\t'<<g.reason<<'\n';
}
void WriteCentromereSupportedIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& s:r.supported_intervals)o<<s.region_id<<'\t'<<s.rb<<'\t'<<s.re<<'\t'<<s.qb<<'\t'<<s.qe<<'\t'<<Sign(s.strand)<<'\n';
}
void WriteCentromereOrientationWindows(std::ostream& o,const CentromereBlockResult& r) {
  o<<"window_id\tquery_begin\tquery_end\tforward_count\treverse_count\tforward_length\treverse_length\tlabel\n";
  for(const auto& w:r.orientation_windows)o<<w.id<<'\t'<<w.begin<<'\t'<<w.end<<'\t'<<w.forward<<'\t'<<w.reverse
    <<'\t'<<w.forward_length<<'\t'<<w.reverse_length<<'\t'<<w.label<<'\n';
}
void WriteCentromereOrientationIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<(x.inferred?"count_inferred_gap":"root_exact_support")<<'\t'<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationSources(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\troot_candidate_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_sources) {
    const auto& c=r.root_candidates.at(x.root_id-1);const auto& a=c.anchor;
    o<<x.interval_id<<'\t'<<x.root_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
      <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(c.strand)<<'\n';
  }
}
void WriteCentromereOrientationGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason";
  if(r.orientation_connected_gaps)o<<"\tleft_block_id\tright_block_id";
  o<<'\n';
  for(const auto& g:r.orientation_gaps) {
    o<<g.begin<<'\t'<<g.end<<'\t'<<Sign(g.strand)<<'\t'<<g.window_id<<'\t'<<g.accepted<<'\t'<<g.reason;
    if(r.orientation_connected_gaps)o<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id;
    o<<'\n';
  }
}
void WriteCentromereOrientationUnknown(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_unknown)o<<x.begin<<'\t'<<x.end<<"\t.\n";
}
void WriteCentromereOrientationBlockAudit(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\tgap_begin\tgap_end\treason\n";
  for(const auto& x:r.orientation_block_audits)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.gap_begin<<'\t'<<x.gap_end<<"\topposite_window_in_old_gap\n";
}
void WriteCentromereOrientationKnownIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_known_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\n';
}
void WriteCentromereOrientationLongGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_long_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<"\tlong_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationLongGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_long_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<Sign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationSparsePivots(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\tparent_begin\tparent_end\treference_status\n";
  for(const auto& x:r.orientation_sparse_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<(x.inferred?"sparse_inferred_gap":"sparse_root_exact")<<'\t'<<x.parent_begin<<'\t'<<x.parent_end
    <<"\tunresolved_copy\n";
}
void WriteCentromereOrientationSparseGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\tparent_begin\tparent_end\taccepted\treason\n";
  for(const auto& g:r.orientation_sparse_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<Sign(g.strand)
    <<'\t'<<g.parent_begin<<'\t'<<g.parent_end<<'\t'<<g.accepted<<'\t'<<g.reason<<'\n';
}
void WriteCentromereOrientationSparseSources(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\troot_candidate_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_sparse_sources) {
    const auto& c=r.root_candidates.at(x.root_id-1);const auto& a=c.anchor;
    o<<x.interval_id<<'\t'<<x.root_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
      <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(c.strand)<<'\n';
  }
}
void WriteCentromereOrientationShortGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_short_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<"\tshort_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationShortGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_short_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<Sign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationMixedGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_mixed_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<"\tmixed_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationMixedGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_mixed_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<Sign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationSparseGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\tparent_begin\tparent_end\tedge_begin\tedge_end\tside\tsource_evidence_id\troot_begin\troot_end\treason\treference_status";
  if(r.prefer_stronger_sparse_gap_end)o<<"\tmeeting_compared\taccepted_root_support\tcompeting_root_support\tscan_begin\tscan_end\tscan_reason";
  o<<'\n';
  for(const auto& x:r.orientation_sparse_gap_intervals) {
    o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<"\tsparse_gap_endpoint\t"<<x.parent_begin<<'\t'<<x.parent_end<<'\t'<<x.edge_begin<<'\t'<<x.edge_end
    <<'\t'<<(x.from_left?"left":"right")<<'\t'<<x.source_interval_id<<'\t'<<x.root_begin<<'\t'<<x.root_end
    <<'\t'<<x.reason<<"\tunresolved_copy";
    if(r.prefer_stronger_sparse_gap_end)o<<'\t'<<x.meeting_compared<<'\t'<<x.accepted_root_support
      <<'\t'<<x.competing_root_support<<'\t'<<x.scan_begin<<'\t'<<x.scan_end<<'\t'<<x.scan_reason;
    o<<'\n';
  }
}
void WriteCentromereOrientationSparseGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\tparent_begin\tparent_end\tedge_begin\tedge_end\tside\tsource_evidence_id\troot_begin\troot_end\taccepted\treason";
  if(r.prefer_stronger_sparse_gap_end)o<<"\tmeeting_compared\taccepted_root_support\tcompeting_root_support\tscan_begin\tscan_end\tscan_reason";
  o<<'\n';
  for(const auto& x:r.orientation_sparse_gap_gaps) {
    o<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.parent_begin<<'\t'<<x.parent_end<<'\t'<<x.edge_begin<<'\t'<<x.edge_end
    <<'\t'<<(x.from_left?"left":"right")<<'\t'<<x.source_interval_id<<'\t'<<x.root_begin<<'\t'<<x.root_end
    <<'\t'<<x.accepted<<'\t'<<x.reason;
    if(r.prefer_stronger_sparse_gap_end)o<<'\t'<<x.meeting_compared<<'\t'<<x.accepted_root_support
      <<'\t'<<x.competing_root_support<<'\t'<<x.scan_begin<<'\t'<<x.scan_end<<'\t'<<x.scan_reason;
    o<<'\n';
  }
}
void WriteCentromereOrientationUnifiedPivots(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tsource_kind\tsource_id\tblock_id\tstrand\tfrom_left\treference_begin\treference_end\tquery_begin\tquery_end\n";
  for(const auto& x:r.orientation_unified_pivots)o<<x.id<<'\t'
    <<(x.source_kind==0?"planning_anchor":x.source_kind==1?"root_exact":"sparse_root_exact")
    <<'\t'<<x.source_id<<'\t'<<x.block_id<<'\t'<<Sign(x.strand)<<'\t'<<x.from_left
    <<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\n';
}
void WriteCentromereOrientationUnifiedScans(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tpivot_id\tstrand\tquery_begin\tquery_end\treason\n";
  for(const auto& x:r.orientation_unified_proposals)o<<x.id<<'\t'<<x.pivot_id<<'\t'<<Sign(x.strand)
    <<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.reason<<'\n';
}
void WriteCentromereOrientationUnifiedEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tproposal_id\n";
  for(const auto& x:r.orientation_unified_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<Sign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereOrientationLocalReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\n";
  for(const auto& x:r.orientation_local_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.left_pivot<<'\t'<<x.right_pivot<<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end
    <<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereOrientationLocalGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  std::uint64_t id=0;
  for(const auto& x:r.orientation_local_reviews)if(x.accepted)o<<++id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<Sign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationIsolatedReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\n";
  for(const auto& x:r.orientation_isolated_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.left_pivot<<'\t'<<x.right_pivot<<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end
    <<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereOrientationOppositeReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\tneighbor_id\tverdict\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_opposite_reviews)o<<x.review_id<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin
    <<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\t'<<x.neighbor_id<<'\t'<<x.verdict<<'\n';
}
void WriteCentromereOrientationIsolatedGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  std::uint64_t id=0;
  for(const auto& x:r.orientation_isolated_reviews)if(x.accepted)o<<++id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<Sign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationCollinearReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\tcompetitive\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_collinear_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.left_pivot<<'\t'<<x.right_pivot<<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end
    <<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.competitive<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\n';
}
void WriteCentromereOrientationCollinearCounters(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\tneighbor_id\tverdict\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_collinear_counters)o<<x.review_id<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin
    <<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\t'<<x.neighbor_id<<'\t'<<x.verdict<<'\n';
}
void WriteCentromereOrientationCollinearWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_collinear_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id
    <<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereOrientationCollinearGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  std::uint64_t id=0;
  for(const auto& x:r.orientation_collinear_reviews)if(x.accepted)o<<++id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<Sign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationAlternateReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\tcompetitive\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_alternate_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.left_pivot<<'\t'<<x.right_pivot<<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end
    <<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.competitive<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\n';
}
void WriteCentromereOrientationAlternateCounters(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\tneighbor_id\tverdict\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_alternate_counters)o<<x.review_id<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin
    <<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\t'<<x.neighbor_id<<'\t'<<x.verdict<<'\n';
}
void WriteCentromereOrientationAlternateWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_alternate_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id
    <<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereOrientationAlternateGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  std::uint64_t id=0;
  for(const auto& x:r.orientation_alternate_reviews)if(x.accepted)o<<++id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<Sign(x.strand)<<'\t'<<x.id<<'\n';
}
} // namespace ramag

namespace ramag {
void WriteCentromereRegionReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_node\tright_node\tcorridor_begin\tcorridor_end\taccepted\treason\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\tsearch_id\n"<<std::setprecision(17);
  for(const auto& x:r.region_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\t'<<x.left_pivot<<'\t'<<x.right_pivot
    <<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end<<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\t'<<x.search_id<<'\n';
}
void WriteCentromereRegionNodes(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\troot_id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\tneutral\n"<<std::setprecision(17);
  for(const auto& x:r.region_nodes){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.root_id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\t'<<x.neutral<<'\n';}
}
void WriteCentromereRegionWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnode_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.region_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereRegionGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  for(const auto& x:r.region_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereRegionSearches(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\treview_id\tleft_node\tright_node\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tinput_bases\tcandidate_count\tsearched\treason\n";
  for(const auto& x:r.region_searches)o<<x.id<<'\t'<<x.review_id<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.reference_begin<<'\t'<<x.reference_end
    <<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<Sign(x.strand)<<'\t'<<x.input_bases<<'\t'<<x.candidate_count<<'\t'<<x.searched<<'\t'<<x.reason<<'\n';
}
} // namespace ramag
namespace ramag {
void WriteCentromereRegionLocalCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\n"<<std::setprecision(17);
  for(const auto& x:r.region_local_candidates){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\n';}
}
} // namespace ramag

namespace ramag {
void WriteCentromereUnifiedReviewReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_node\tright_node\tcorridor_begin\tcorridor_end\taccepted\treason\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\tsearch_id\n"<<std::setprecision(17);
  for(const auto& x:r.unified_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\t'<<x.left_pivot<<'\t'<<x.right_pivot
    <<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end<<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\t'<<x.search_id<<'\n';
}
void WriteCentromereUnifiedReviewNodes(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\troot_id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\tneutral\n"<<std::setprecision(17);
  for(const auto& x:r.unified_nodes){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.root_id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\t'<<x.neutral<<'\n';}
}
void WriteCentromereUnifiedReviewWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnode_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.unified_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereUnifiedReviewIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  for(const auto& x:r.unified_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereUnifiedReviewSearches(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\treview_id\tleft_node\tright_node\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tinput_bases\tcandidate_count\tsearched\treason\n";
  for(const auto& x:r.unified_searches)o<<x.id<<'\t'<<x.review_id<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.reference_begin<<'\t'<<x.reference_end
    <<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<Sign(x.strand)<<'\t'<<x.input_bases<<'\t'<<x.candidate_count<<'\t'<<x.searched<<'\t'<<x.reason<<'\n';
}
} // namespace ramag
namespace ramag {
void WriteCentromereUnifiedReviewLocalCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\n"<<std::setprecision(17);
  for(const auto& x:r.unified_local_candidates){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<Sign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\n';}
}
void WriteCentromereUnifiedReviewGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\troot_review_id\tchosen_review_id\thypothesis_count\taccepted\treason\n";
  for(const auto& x:r.unified_gap_decisions)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)
    <<'\t'<<x.root_review_id<<'\t'<<x.chosen_review_id<<'\t'<<x.hypothesis_count<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereUnifiedReviewHypotheses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tgap_id\tleft_node\tright_node\tsearch_id\treview_id\taccepted\treason\n";
  for(const auto& x:r.unified_hypotheses)o<<x.id<<'\t'<<x.gap_id<<'\t'<<x.left_node<<'\t'<<x.right_node
    <<'\t'<<x.search_id<<'\t'<<x.review_id<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
} // namespace ramag

namespace ramag {
void WriteCentromereKnownCounterchainReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tblock_id\tquery_begin\tquery_end\tgap_begin\tgap_end\tstrand\tleft_node\tright_node\topposite_support\topposite_score\topposite_nodes\town_support\town_score\town_nodes\tconnection_checks\tretracted\treason\n"<<std::setprecision(17);
  for(const auto& x:r.known_counterchain_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.gap_begin<<'\t'<<x.gap_end<<'\t'<<Sign(x.strand)<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.opposite_support<<'\t'<<x.opposite_score<<'\t'<<x.opposite_nodes<<'\t'<<x.own_support<<'\t'<<x.own_score<<'\t'<<x.own_nodes<<'\t'<<x.connection_checks<<'\t'<<x.retracted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereKnownCounterchainWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.known_counterchain_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereKnownCounterchainIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.known_counterchain_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\n';
}
} // namespace ramag

namespace ramag {
void WriteCentromereLocalDirectionRepairReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tblock_id\tquery_begin\tquery_end\tstrand\tinherited_strand\tsupport\topposite_support\tscore\topposite_score\tnodes\topposite_nodes\taccepted\treason\n"<<std::setprecision(17);
  for(const auto& x:r.local_direction_repair_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\t'<<Sign(x.inherited_strand)<<'\t'<<x.support<<'\t'<<x.opposite_support<<'\t'<<x.score<<'\t'<<x.opposite_score<<'\t'<<x.nodes<<'\t'<<x.opposite_nodes<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereLocalDirectionRepairWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.local_direction_repair_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereLocalDirectionRepairIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.local_direction_repair_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<Sign(x.strand)<<'\n';
}
void WriteCentromereLocalDirectionRepairChanges(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tquery_begin\tquery_end\tbefore_strand\tafter_strand\treason\n";
  for(const auto& x:r.local_direction_repair_changes)o<<x.review_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.before_strand<<'\t'<<x.after_strand<<'\t'<<x.reason<<'\n';
}
} // namespace ramag

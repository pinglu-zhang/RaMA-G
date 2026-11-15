#include "ramag/centromere_blocks.hpp"
#include "ramag/internal/centromere/refinement.hpp"
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include "ramag/centromere_direction.hpp"
#include "ramag/internal/common/centromere.hpp"
#include "ramag/centromere_stitch_formats.hpp"
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "ramag/internal/centromere/gap_alignment.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <new>
#include <tuple>
#include <unordered_map>
namespace ramag {
void WriteCentromereDirectionIslandContextReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tblock_id\tchecks\tcomponents\toriginal_begin\toriginal_end\tcore_reference_begin\tcore_reference_end\tcore_query_begin\tcore_query_end\treference_explained\tquery_explained\treference_unsupported\tquery_unsupported\told_score\tscore\taccepted\treason\n";
  for(const auto& x:r.direction_island_context_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.checks<<'\t'<<x.components<<'\t'
    <<x.original_begin<<'\t'<<x.original_end<<'\t'<<x.core_reference_begin<<'\t'<<x.core_reference_end<<'\t'
    <<x.core_query_begin<<'\t'<<x.core_query_end<<'\t'<<x.reference_explained<<'\t'<<x.query_explained<<'\t'
    <<x.reference_unsupported<<'\t'<<x.query_unsupported<<'\t'<<x.old_score<<'\t'<<x.score<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
  if(r.direction_island_context_budget_exhausted)o<<"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tremaining_batch_checks_exhausted\n";
}
void WriteCentromereDirectionPlanReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tblock_id\tleft_index\tsource_id\tchecks\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tweight\taccepted\treason\n";
  for(const auto& x:r.direction_plan_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.left_index<<'\t'<<x.source_id<<'\t'<<x.checks<<'\t'
    <<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<(x.strand==Strand::Forward?'+':'-')
    <<'\t'<<x.weight<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
  if(r.direction_plan_review_budget_exhausted)o<<"0\t0\t0\t0\t"<<r.direction_plan_review_checks<<"\t0\t0\t0\t0\t+\t0\t0\tbatch_budget_exhausted\n";
}
void WriteCentromereDirectionLocalProtectionReviews(std::ostream& o,const CentromereBlockResult& result) {
  o<<"id\tblock_id\tbegin\tend\tstrand\told_matched\tcounter_matched\tcandidate_same\tcandidate_opposite\tcandidate_unknown\trestored\tprotected_interval\tregion_changed\treason\n";
  for(const auto& a:result.direction_local_protection_reviews)o<<a.id<<'\t'<<a.block_id<<'\t'<<a.begin<<'\t'<<a.end<<'\t'
    <<(a.strand==Strand::Forward?'+':'-')<<'\t'<<a.old_matched<<'\t'<<a.counter_matched<<'\t'<<a.candidate_same<<'\t'
    <<a.candidate_opposite<<'\t'<<a.candidate_unknown<<'\t'<<a.restored<<'\t'<<a.protected_interval<<'\t'<<a.region_changed<<'\t'<<a.reason<<'\n';
}
void WriteCentromereDirectionIslandAudits(std::ostream& o,const CentromereBlockResult& result) {
  o<<"id\tblock_id\tleft_source_id\tright_source_id\tchecks\tbegin\tend\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tbackground\tweight\treference_covered\tquery_covered\treference_unsupported\tquery_unsupported\tscore\teligible\tretracted\treason\tsource_role\n";
  for(const auto& a:result.direction_island_audits)o<<a.id<<'\t'<<a.block_id<<'\t'<<a.left_source_id<<'\t'<<a.right_source_id<<'\t'<<a.checks<<'\t'
    <<a.begin<<'\t'<<a.end<<'\t'<<a.reference_begin<<'\t'<<a.reference_end<<'\t'<<a.query_begin<<'\t'<<a.query_end<<'\t'
    <<(a.strand==Strand::Forward?'+':'-')<<'\t'<<(a.background==Strand::Forward?'+':'-')<<'\t'<<a.weight<<'\t'
    <<a.reference_covered<<'\t'<<a.query_covered<<'\t'<<a.reference_unsupported<<'\t'<<a.query_unsupported<<'\t'<<a.score<<'\t'
    <<a.eligible<<'\t'<<a.retracted<<'\t'<<a.reason<<'\t'<<a.source_role<<'\n';
}
void WriteCentromereOverlapVariants(std::ostream& o,const CentromereBlockResult& r) {
  o<<"variant_id\tstrand\toriginal_reference_begin\toriginal_query_begin\toriginal_length\treference_begin\tquery_begin\tlength\ttrim_front\ttrim_back\n";
  std::uint64_t id=0;
  for(const auto& v:r.overlap_variants)o<<++id<<'\t'<<(v.strand==Strand::Forward?'+':'-')<<'\t'
    <<v.original.reference_begin<<'\t'<<v.original.query_begin<<'\t'<<v.original.length<<'\t'
    <<v.variant.reference_begin<<'\t'<<v.variant.query_begin<<'\t'<<v.variant.length<<'\t'
    <<v.trim_front<<'\t'<<v.trim_back<<'\n';
}
void WriteCentromereOverlapSummary(std::ostream& o,const CentromereBlockResult& r) {
  o<<"pairs_examined\tvariants_added\tbudget_exhausted\tnonuniform_nodes_skipped\tlocal_support_rejections\n"<<r.overlap_pairs_examined<<'\t'
    <<r.overlap_variants_added<<'\t'<<r.overlap_budget_exhausted<<'\t'<<r.overlap_nonuniform_nodes_skipped<<'\t'<<r.overlap_local_support_rejections<<'\n';
}
}
namespace ramag {
namespace {
char BlockSign(Strand strand) {return strand==Strand::Forward ? '+' : '-';}
void Header(std::ostream& o,bool joint) {o<<"block_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tskeleton\ttentative\tinversion_candidate\tforward_support\treverse_support\tstatus\treason";if(joint)o<<"\tdirection_known";o<<'\n';}
void Row(std::ostream& o,const CentromereBlock& b,bool joint) {
  o<<b.id<<'\t'<<b.reference_begin<<'\t'<<b.reference_end<<'\t'<<b.query_begin<<'\t'<<b.query_end<<'\t'<<(b.direction_known ? BlockSign(b.strand) : '.')<<'\t'
   <<b.skeleton<<'\t'<<b.tentative<<'\t'<<b.inversion_candidate<<'\t'<<b.forward_support<<'\t'<<b.reverse_support<<'\t'<<b.status<<'\t'<<b.reason;if(joint)o<<'\t'<<b.direction_known;o<<'\n';
}
void AnchorHeader(std::ostream& o) {o<<"anchor_id\tblock_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tmatch_length\tweight\n";o<<std::setprecision(std::numeric_limits<double>::max_digits10);}
void AnchorRow(std::ostream& o,const CentromereAnchor& a,Strand strand,std::uint64_t id,std::uint64_t block) {
  o<<id<<'\t'<<block<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\n';
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
     <<'\t'<<BlockSign(c.strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\t'<<c.reference_count<<'\t'<<c.query_count<<'\t'<<c.rare_level<<'\n';
  }
}
void WriteCentromereDirectionRegions(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\texact_support\tsecondary_score\tcompeting_mappings\tblock_count\tleft_boundary_width\tright_boundary_width\tpaired_columns\tinsertions\tdeletions\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(const auto& x:r.regions)o<<x.id<<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<BlockSign(x.strand)<<'\t'
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
      <<'\t'<<BlockSign(c.strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\t'<<c.reference_count<<'\t'<<c.query_count<<'\t'<<c.rare_level<<'\t'
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
      <<'\t'<<BlockSign(n.strand)<<'\t'<<a.weight<<'\t'<<n.predecessor<<'\t'<<n.first<<'\t'<<n.count<<'\t'<<n.support<<'\t'<<n.weight<<'\t'<<n.neutral<<'\n';
  }
}
void WriteCentromereCountCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"candidate_id\tfirst\tendpoint\tcount\tforward_count\treverse_count\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\texact_support\tweight\teligible\tselected\n";
  o<<std::setprecision(std::numeric_limits<double>::max_digits10);
  for(const auto& c:r.count_candidates)o<<c.id<<'\t'<<c.first<<'\t'<<c.endpoint<<'\t'<<c.count<<'\t'<<c.forward<<'\t'<<c.reverse
    <<'\t'<<c.rb<<'\t'<<c.re<<'\t'<<c.qb<<'\t'<<c.qe<<'\t'<<BlockSign(c.strand)<<'\t'<<c.support<<'\t'<<c.weight<<'\t'<<c.eligible<<'\t'<<c.selected<<'\n';
}
void WriteCentromereCountGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\tleft_node\tright_node\twindow_id\tforward_count\treverse_count\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\taccepted\treason\n";
  for(const auto& g:r.count_gaps)o<<g.region_id<<'\t'<<g.left_node<<'\t'<<g.right_node<<'\t'<<g.window_id<<'\t'<<g.forward<<'\t'<<g.reverse
    <<'\t'<<g.rb<<'\t'<<g.re<<'\t'<<g.qb<<'\t'<<g.qe<<'\t'<<BlockSign(g.strand)<<'\t'<<g.accepted<<'\t'<<g.reason<<'\n';
}
void WriteCentromereSupportedIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& s:r.supported_intervals)o<<s.region_id<<'\t'<<s.rb<<'\t'<<s.re<<'\t'<<s.qb<<'\t'<<s.qe<<'\t'<<BlockSign(s.strand)<<'\n';
}
void WriteCentromereOrientationWindows(std::ostream& o,const CentromereBlockResult& r) {
  o<<"window_id\tquery_begin\tquery_end\tforward_count\treverse_count\tforward_length\treverse_length\tlabel\n";
  for(const auto& w:r.orientation_windows)o<<w.id<<'\t'<<w.begin<<'\t'<<w.end<<'\t'<<w.forward<<'\t'<<w.reverse
    <<'\t'<<w.forward_length<<'\t'<<w.reverse_length<<'\t'<<w.label<<'\n';
}
void WriteCentromereOrientationIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<'\t'<<(x.inferred?"count_inferred_gap":"root_exact_support")<<'\t'<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationSources(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\troot_candidate_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_sources) {
    const auto& c=r.root_candidates.at(x.root_id-1);const auto& a=c.anchor;
    o<<x.interval_id<<'\t'<<x.root_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
      <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(c.strand)<<'\n';
  }
}
void WriteCentromereOrientationGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason";
  if(r.orientation_connected_gaps)o<<"\tleft_block_id\tright_block_id";
  o<<'\n';
  for(const auto& g:r.orientation_gaps) {
    o<<g.begin<<'\t'<<g.end<<'\t'<<BlockSign(g.strand)<<'\t'<<g.window_id<<'\t'<<g.accepted<<'\t'<<g.reason;
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
  for(const auto& x:r.orientation_block_audits)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<'\t'<<x.gap_begin<<'\t'<<x.gap_end<<"\topposite_window_in_old_gap\n";
}
void WriteCentromereOrientationKnownIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_known_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\n';
}
void WriteCentromereOrientationLongGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_long_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<"\tlong_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationLongGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_long_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<BlockSign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationSparsePivots(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\tparent_begin\tparent_end\treference_status\n";
  for(const auto& x:r.orientation_sparse_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<'\t'<<(x.inferred?"sparse_inferred_gap":"sparse_root_exact")<<'\t'<<x.parent_begin<<'\t'<<x.parent_end
    <<"\tunresolved_copy\n";
}
void WriteCentromereOrientationSparseGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\tparent_begin\tparent_end\taccepted\treason\n";
  for(const auto& g:r.orientation_sparse_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<BlockSign(g.strand)
    <<'\t'<<g.parent_begin<<'\t'<<g.parent_end<<'\t'<<g.accepted<<'\t'<<g.reason<<'\n';
}
void WriteCentromereOrientationSparseSources(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\troot_candidate_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.orientation_sparse_sources) {
    const auto& c=r.root_candidates.at(x.root_id-1);const auto& a=c.anchor;
    o<<x.interval_id<<'\t'<<x.root_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
      <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(c.strand)<<'\n';
  }
}
void WriteCentromereOrientationShortGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_short_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<"\tshort_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationShortGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_short_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<BlockSign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationMixedGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\twindow_id\treference_status\n";
  for(const auto& x:r.orientation_mixed_gap_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<"\tmixed_gap_endpoint\t"<<x.window_id<<"\tunresolved_copy\n";
}
void WriteCentromereOrientationMixedGapGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"query_begin\tquery_end\tstrand\twindow_id\taccepted\treason\tleft_block_id\tright_block_id\n";
  for(const auto& g:r.orientation_mixed_gap_gaps)o<<g.begin<<'\t'<<g.end<<'\t'<<BlockSign(g.strand)<<'\t'<<g.window_id
    <<'\t'<<g.accepted<<'\t'<<g.reason<<'\t'<<g.left_block_id<<'\t'<<g.right_block_id<<'\n';
}
void WriteCentromereOrientationSparseGapEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"evidence_id\tquery_begin\tquery_end\tstrand\tsource\tparent_begin\tparent_end\tedge_begin\tedge_end\tside\tsource_evidence_id\troot_begin\troot_end\treason\treference_status";
  if(r.prefer_stronger_sparse_gap_end)o<<"\tmeeting_compared\taccepted_root_support\tcompeting_root_support\tscan_begin\tscan_end\tscan_reason";
  o<<'\n';
  for(const auto& x:r.orientation_sparse_gap_intervals) {
    o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
    o<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
    <<'\t'<<x.source_id<<'\t'<<x.block_id<<'\t'<<BlockSign(x.strand)<<'\t'<<x.from_left
    <<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\n';
}
void WriteCentromereOrientationUnifiedScans(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tpivot_id\tstrand\tquery_begin\tquery_end\treason\n";
  for(const auto& x:r.orientation_unified_proposals)o<<x.id<<'\t'<<x.pivot_id<<'\t'<<BlockSign(x.strand)
    <<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.reason<<'\n';
}
void WriteCentromereOrientationUnifiedEnds(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tproposal_id\n";
  for(const auto& x:r.orientation_unified_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<BlockSign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereOrientationLocalReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\n";
  for(const auto& x:r.orientation_local_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
    <<'\t'<<x.left_pivot<<'\t'<<x.right_pivot<<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end
    <<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereOrientationLocalGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  std::uint64_t id=0;
  for(const auto& x:r.orientation_local_reviews)if(x.accepted)o<<++id<<'\t'<<x.begin<<'\t'<<x.end
    <<'\t'<<BlockSign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationIsolatedReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\n";
  for(const auto& x:r.orientation_isolated_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
    <<'\t'<<BlockSign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationCollinearReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\tcompetitive\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_collinear_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
    <<'\t'<<BlockSign(x.strand)<<'\t'<<x.id<<'\n';
}
void WriteCentromereOrientationAlternateReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_pivot_id\tright_pivot_id\tcorridor_begin\tcorridor_end\taccepted\treason\tcompetitive\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\n";
  o<<std::setprecision(17);
  for(const auto& x:r.orientation_alternate_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
    <<'\t'<<BlockSign(x.strand)<<'\t'<<x.id<<'\n';
}
} // namespace ramag

namespace ramag {
void WriteCentromereRegionReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_node\tright_node\tcorridor_begin\tcorridor_end\taccepted\treason\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\tsearch_id\n"<<std::setprecision(17);
  for(const auto& x:r.region_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.left_pivot<<'\t'<<x.right_pivot
    <<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end<<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\t'<<x.search_id<<'\n';
}
void WriteCentromereRegionNodes(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\troot_id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\tneutral\n"<<std::setprecision(17);
  for(const auto& x:r.region_nodes){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.root_id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\t'<<x.neutral<<'\n';}
}
void WriteCentromereRegionWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnode_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.region_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereRegionGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  for(const auto& x:r.region_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereRegionSearches(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\treview_id\tleft_node\tright_node\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tinput_bases\tcandidate_count\tsearched\treason\n";
  for(const auto& x:r.region_searches)o<<x.id<<'\t'<<x.review_id<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.reference_begin<<'\t'<<x.reference_end
    <<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.input_bases<<'\t'<<x.candidate_count<<'\t'<<x.searched<<'\t'<<x.reason<<'\n';
}
} // namespace ramag
namespace ramag {
void WriteCentromereRegionLocalCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\n"<<std::setprecision(17);
  for(const auto& x:r.region_local_candidates){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\n';}
}
} // namespace ramag

namespace ramag {
void WriteCentromereUnifiedReviewReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\tleft_node\tright_node\tcorridor_begin\tcorridor_end\taccepted\treason\tcandidate_nodes\tconnection_checks\tchain_nodes\tsupport\tchain_score\tsearch_id\n"<<std::setprecision(17);
  for(const auto& x:r.unified_reviews)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.left_pivot<<'\t'<<x.right_pivot
    <<'\t'<<x.corridor_begin<<'\t'<<x.corridor_end<<'\t'<<x.accepted<<'\t'<<x.reason<<'\t'<<x.candidate_nodes<<'\t'<<x.connection_checks
    <<'\t'<<x.chain_nodes<<'\t'<<x.support<<'\t'<<x.chain_score<<'\t'<<x.search_id<<'\n';
}
void WriteCentromereUnifiedReviewNodes(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\troot_id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\tneutral\n"<<std::setprecision(17);
  for(const auto& x:r.unified_nodes){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.root_id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\t'<<x.neutral<<'\n';}
}
void WriteCentromereUnifiedReviewWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnode_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.unified_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereUnifiedReviewIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\treview_id\n";
  for(const auto& x:r.unified_intervals)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.scan_id<<'\n';
}
void WriteCentromereUnifiedReviewSearches(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\treview_id\tleft_node\tright_node\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tinput_bases\tcandidate_count\tsearched\treason\n";
  for(const auto& x:r.unified_searches)o<<x.id<<'\t'<<x.review_id<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.reference_begin<<'\t'<<x.reference_end
    <<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.input_bases<<'\t'<<x.candidate_count<<'\t'<<x.searched<<'\t'<<x.reason<<'\n';
}
} // namespace ramag
namespace ramag {
void WriteCentromereUnifiedReviewLocalCandidates(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tsearch_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\treference_count\tquery_count\trare_level\n"<<std::setprecision(17);
  for(const auto& x:r.unified_local_candidates){const auto& a=x.source.anchor;o<<x.id<<'\t'<<x.search_id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length
    <<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'<<BlockSign(x.source.strand)<<'\t'<<a.length<<'\t'<<a.weight
    <<'\t'<<x.source.reference_count<<'\t'<<x.source.query_count<<'\t'<<x.source.rare_level<<'\n';}
}
void WriteCentromereUnifiedReviewGaps(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tquery_begin\tquery_end\tstrand\troot_review_id\tchosen_review_id\thypothesis_count\taccepted\treason\n";
  for(const auto& x:r.unified_gap_decisions)o<<x.id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)
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
  for(const auto& x:r.known_counterchain_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.gap_begin<<'\t'<<x.gap_end<<'\t'<<BlockSign(x.strand)<<'\t'<<x.left_node<<'\t'<<x.right_node<<'\t'<<x.opposite_support<<'\t'<<x.opposite_score<<'\t'<<x.opposite_nodes<<'\t'<<x.own_support<<'\t'<<x.own_score<<'\t'<<x.own_nodes<<'\t'<<x.connection_checks<<'\t'<<x.retracted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereKnownCounterchainWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.known_counterchain_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereKnownCounterchainIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.known_counterchain_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\n';
}
} // namespace ramag

namespace ramag {
void WriteCentromereLocalDirectionRepairReviews(std::ostream& o,const CentromereBlockResult& r) {
  o<<"id\tblock_id\tquery_begin\tquery_end\tstrand\tinherited_strand\tsupport\topposite_support\tscore\topposite_score\tnodes\topposite_nodes\taccepted\treason\n"<<std::setprecision(17);
  for(const auto& x:r.local_direction_repair_reviews)o<<x.id<<'\t'<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\t'<<BlockSign(x.inherited_strand)<<'\t'<<x.support<<'\t'<<x.opposite_support<<'\t'<<x.score<<'\t'<<x.opposite_score<<'\t'<<x.nodes<<'\t'<<x.opposite_nodes<<'\t'<<x.accepted<<'\t'<<x.reason<<'\n';
}
void WriteCentromereLocalDirectionRepairWitnesses(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tordinal\tnormalized_id\treference_begin\tquery_begin\tlength\tweight\n"<<std::setprecision(17);
  for(const auto& x:r.local_direction_repair_witnesses)o<<x.review_id<<'\t'<<x.ordinal<<'\t'<<x.normalized_id<<'\t'<<x.reference_begin<<'\t'<<x.query_begin<<'\t'<<x.length<<'\t'<<x.weight<<'\n';
}
void WriteCentromereLocalDirectionRepairIntervals(std::ostream& o,const CentromereBlockResult& r) {
  o<<"block_id\tquery_begin\tquery_end\tstrand\n";
  for(const auto& x:r.local_direction_repair_intervals)o<<x.block_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<BlockSign(x.strand)<<'\n';
}
void WriteCentromereLocalDirectionRepairChanges(std::ostream& o,const CentromereBlockResult& r) {
  o<<"review_id\tquery_begin\tquery_end\tbefore_strand\tafter_strand\treason\n";
  for(const auto& x:r.local_direction_repair_changes)o<<x.review_id<<'\t'<<x.begin<<'\t'<<x.end<<'\t'<<x.before_strand<<'\t'<<x.after_strand<<'\t'<<x.reason<<'\n';
}
} // namespace ramag


namespace ramag {
namespace {
char DirectionSign(Strand s) {return s==Strand::Forward ? '+' : '-';}
void Box(std::ostream& o,const DirectionBox& b) {
  o<<b.reference_begin<<'\t'<<b.reference_end<<'\t'<<b.query_begin<<'\t'<<b.query_end;
}
void Precision(std::ostream& o) {o<<std::setprecision(std::numeric_limits<double>::max_digits10);}
void List(std::ostream& o,const std::vector<std::uint64_t>& ids) {
  for(std::size_t i=0;i<ids.size();++i) {if(i)o<<',';o<<ids[i];}
}
}
void WriteDirectionAnchors(std::ostream& o,const CentromereDirectionResult& r) {
  Precision(o);
  o<<"anchor_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tsource_weight\tfull_reference_count\tfull_query_count\ttask_id\tdepth\trare_level\treference_interval_count\tquery_interval_count\tfirst_hit_level\n";
  for(const auto& a:r.anchors) for(const auto& p:a.origins)
    o<<a.id<<'\t'<<a.reference_begin<<'\t'<<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'
     <<DirectionSign(a.strand)<<'\t'<<a.length<<'\t'<<a.source_weight<<'\t'<<a.full_reference_count<<'\t'<<a.full_query_count<<'\t'
     <<p.task_id<<'\t'<<p.depth<<'\t'<<p.rare_level<<'\t'<<p.reference_interval_count<<'\t'<<p.query_interval_count<<'\t'<<p.first_hit_level<<'\n';
}
void WriteDirectionRegions(std::ostream& o,const CentromereDirectionResult& r) {
  Precision(o);
  o<<"region_id\ttask_id\tdepth\treference_begin\treference_end\tquery_begin\tquery_end\tscale\tforward_support\treverse_support\tadvantage\tdirection\tforward_longest\treverse_longest\tforward_chain_nodes\treverse_chain_nodes\n";
  for(const auto& a:r.regions) {
    o<<a.id<<'\t'<<a.task_id<<'\t'<<a.depth<<'\t';Box(o,a.box);
    const auto f=a.chains[0].support,b=a.chains[1].support;
    o<<'\t'<<a.scale<<'\t'<<f<<'\t'<<b<<'\t'<<a.advantage<<'\t'
     <<(!f&&!b ? "insufficient" : f==b ? "tied" : f>b ? "forward_preferred" : "reverse_preferred");
    std::array<Length,2> longest{};std::array<std::size_t,2> count{};
    for(std::size_t d=0;d<2;++d)
      for(auto i=a.chains[d].best;i!=kNoDirectionNode;i=a.chains[d].nodes[i].predecessor) {
        longest[d]=std::max(longest[d],a.chains[d].nodes[i].length);++count[d];
      }
    o<<'\t'<<longest[0]<<'\t'<<longest[1]<<'\t'<<count[0]<<'\t'<<count[1]<<'\n';
  }
}
void WriteDirectionChains(std::ostream& o,const CentromereDirectionResult& r) {
  Precision(o);
  o<<"region_id\tstrand\tnode_id\tpredecessor\tbest_endpoint\treference_begin\treference_end\tquery_begin\tquery_end\tlength\tpath_support\tsource_score\tsource_anchor_ids\n";
  for(const auto& region:r.regions) for(std::size_t d=0;d<2;++d) {
    const auto& c=region.chains[d];
    for(std::size_t i=0;i<c.nodes.size();++i) {
      const auto& n=c.nodes[i];
      const auto qb=d==0 ? n.oriented_query_begin : r.query_length-n.oriented_query_begin-n.length;
      o<<region.id<<'\t'<<(d==0 ? '+' : '-')<<'\t'<<i<<'\t';
      if(n.predecessor==kNoDirectionNode)o<<'.';else o<<n.predecessor;
      o<<'\t'<<(i==c.best)<<'\t'<<n.reference_begin<<'\t'<<n.reference_begin+n.length<<'\t'<<qb<<'\t'<<qb+n.length
       <<'\t'<<n.length<<'\t'<<n.support<<'\t'<<n.source_score<<'\t';List(o,n.source_ids);o<<'\n';
    }
  }
}
void WriteDirectionTasks(std::ostream& o,const CentromereDirectionResult& r) {
  o<<"task_id\tparent_id\tdepth\treference_begin\treference_end\tquery_begin\tquery_end\tstop_reason\n";
  for(const auto& a:r.tasks) {o<<a.id<<'\t'<<a.parent_id<<'\t'<<a.depth<<'\t';Box(o,a.box);o<<'\t'<<a.stop_reason<<'\n';}
}
void WriteInversionEvidence(std::ostream& o,const CentromereDirectionResult& r) {
  o<<"candidate_id\tregion_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tbackground\tclassification\treason\tsupport\topposite_support\tcontext_background_support\tcompeting\tleft_flank\tright_flank\tref_left_low\tref_left_high\tref_right_low\tref_right_high\tquery_left_low\tquery_left_high\tquery_right_low\tquery_right_high\tanchor_ids\tindependent_segments\tlongest_match\tenclosing_opposite_candidates\n";
  for(const auto& a:r.candidates) {
    o<<a.id<<'\t'<<a.region_id<<'\t';Box(o,a.block);
    o<<'\t'<<DirectionSign(a.strand)<<'\t'<<DirectionSign(a.background)<<'\t'<<a.classification<<'\t'<<a.reason<<'\t'
     <<a.support<<'\t'<<a.opposite_support<<'\t'<<a.context_background_support<<'\t'<<a.competing<<'\t'<<a.left_flank<<'\t'<<a.right_flank;
    auto bounds=[&](bool known,const auto& interval){if(known)o<<'\t'<<interval[0]<<'\t'<<interval[1];else o<<"\t.\t.";};
    bounds(a.left_bounded,a.reference_left);bounds(a.right_bounded,a.reference_right);
    bounds(a.left_bounded,a.query_left);bounds(a.right_bounded,a.query_right);
    o<<'\t';List(o,a.anchors);o<<'\t'<<a.independent_segments<<'\t'<<a.longest_match<<'\t';List(o,a.enclosing_opposite_candidates);o<<'\n';
  }
}
} // namespace ramag



namespace ramag {
namespace detail {
// Private helper is separately testable without allocating a 256 MiB sequence.
void WriteStitchSamCigarRun(std::ostream& out, Length length, char operation) {
  constexpr Length max_run = (Length{1} << 28U) - 1U;
  while (length) {
    const auto chunk = std::min(length, max_run);
    out << chunk << operation;
    length -= chunk;
  }
}
}  // namespace detail
namespace {
struct Counts {
  std::uint64_t matches{}, paired{}, columns{}, edits{};
};
struct ExportSegment {
  const CentromereStitchSegment* segment{};
  Counts counts;
};

using ::ramag::detail::IsCanonicalUpper;
char Complement(char c) {
  switch (c) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    default: return 'N';
  }
}
char Sign(Strand strand) { return strand == Strand::Forward ? '+' : '-'; }
void CheckName(const std::string& name) {
  if (name.empty() || name == "*") throw AlignmentError("stitch export requires a nonempty sequence name other than '*'");
  for (unsigned char c : name)
    if (c < 33 || c > 126) throw AlignmentError("stitch export sequence names must be printable ASCII without whitespace");
}
void CheckSamNames(const SequenceRecord& reference, const SequenceRecord& query) {
  // SAM QNAME excludes '@'; reference punctuation is restricted by SAM v1.6.
  if (query.name.size() > 254 || query.name.find('@') != std::string::npos)
    throw AlignmentError("stitch SAM query name exceeds 254 characters or contains '@'");
  if (reference.name.front() == '*' || reference.name.front() == '=')
    throw AlignmentError("stitch SAM reference name must not start with '*' or '='");
  constexpr std::string_view punctuation = "!#$%&*+./:;=?@^_|~-";
  for (unsigned char c : reference.name)
    if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) &&
        punctuation.find(static_cast<char>(c)) == std::string_view::npos)
      throw AlignmentError("stitch SAM reference name contains reserved punctuation");
  if (reference.size() == 0 || reference.size() > static_cast<Length>(std::numeric_limits<std::int32_t>::max()))
    throw AlignmentError("stitch SAM reference length is outside the SAM @SQ range");
}
void Add(std::uint64_t& destination, std::uint64_t n) {
  if (n > std::numeric_limits<std::uint64_t>::max() - destination)
    throw AlignmentError("stitch export CIGAR count overflow");
  destination += n;
}
std::vector<ExportSegment> Validate(const SequenceRecord& reference, const SequenceRecord& query,
                                  const CentromereStitchResult& result) {
  if (result.aborted) throw AlignmentError("cannot export an aborted stitch snapshot as alignment results");
  CheckName(reference.name); CheckName(query.name);
  for (const auto* sequence : {&reference, &query})
    for (char c : sequence->bases)
      if (!IsCanonicalUpper(c) && c != 'N') throw AlignmentError("stitch export requires normalized uppercase ACGTN sequence");
  if (result.reference_length != reference.size() || result.query_length != query.size())
    throw AlignmentError("stitch export sequence lengths disagree with stitch result");
  if (result.score_model != CentromereStitchScoreModel::Centrolign3Piece &&
      result.score_model != CentromereStitchScoreModel::Rama2PieceControl)
    throw AlignmentError("invalid stitch export score model");
  std::set<std::uint64_t> ids;
  std::vector<ExportSegment> entries;
  entries.reserve(result.segments.size());
  for (const auto& segment : result.segments) {
    const auto& a = segment.alignment;
    if (!segment.segment_id || !segment.path_id || !ids.insert(segment.segment_id).second ||
        (segment.path_role != "primary" && segment.path_role != "copy"))
      throw AlignmentError("invalid stitch export segment identity or role");
    if (a.reference_id != reference.numeric_id || a.query_id != query.numeric_id ||
        a.reference_begin >= a.reference_end || a.reference_end > reference.size() ||
        a.query_begin >= a.query_end || a.query_end > query.size() || a.cigar.empty() ||
        (a.strand != Strand::Forward && a.strand != Strand::Reverse))
      throw AlignmentError("invalid stitch export segment sequence identity, bounds or strand");
    Position r = a.reference_begin, q = 0;
    Counts counts;
    const Length qspan = a.query_end - a.query_begin;
    for (const auto& op : a.cigar) {
      if (!op.length || (op.operation != '=' && op.operation != 'X' && op.operation != 'I' && op.operation != 'D'))
        throw AlignmentError("stitch export requires nonempty =/X/I/D CIGAR runs");
      if (op.operation != 'I' && op.length > a.reference_end - r)
        throw AlignmentError("stitch export CIGAR overruns reference interval");
      if (op.operation != 'D' && op.length > qspan - q)
        throw AlignmentError("stitch export CIGAR overruns query interval");
      Add(counts.columns, op.length);
      if (op.operation == '=' || op.operation == 'X') {
        Add(counts.paired, op.length);
        for (Length k = 0; k < op.length; ++k) {
          const char rb = reference.bases[r + k];
          const char qb = a.strand == Strand::Forward ? query.bases[a.query_begin + q + k]
              : Complement(query.bases[a.query_end - q - k - 1]);
          if ((IsCanonicalUpper(rb) && rb == qb) != (op.operation == '='))
            throw AlignmentError("stitch export =/X disagrees with sequence characters");
        }
      }
      if (op.operation == '=') Add(counts.matches, op.length);
      else Add(counts.edits, op.length);
      if (op.operation != 'I') r += op.length;
      if (op.operation != 'D') q += op.length;
    }
    if (r != a.reference_end || q != qspan || counts.edits != a.edit_distance)
      throw AlignmentError("stitch export CIGAR consumption or edit distance disagrees with segment");
    entries.push_back({&segment, counts});
  }
  std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
    return a.segment->segment_id < b.segment->segment_id;
  });
  return entries;
}
void Cigar(std::ostream& out, const std::vector<CigarOp>& cigar) {
  for (const auto& op : cigar) out << op.length << op.operation;
}
void IdList(std::ostream& out, const std::vector<std::uint64_t>& ids) {
  for (std::size_t i = 0; i < ids.size(); ++i) { if (i) out << ','; out << ids[i]; }
}
void Tags(std::ostream& out, const ExportSegment& entry, const CentromereStitchResult& result) {
  const auto& s = *entry.segment;
  out << "\tsi:Z:" << s.segment_id << "\tpi:Z:" << s.path_id << "\tpr:Z:" << s.path_role
      << "\tri:Z:" << s.region_id << "\tbi:Z:" << s.branch_id
      << "\tsm:Z:" << CentromereStitchScoreModelName(result.score_model)
      << "\tst:Z:" << (result.complete ? "complete" : "partial") << "\tsc:Z:" << s.alignment.score;
  if (entry.counts.edits <= static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()))
    out << "\tNM:i:" << entry.counts.edits;
  else out << "\tne:Z:" << entry.counts.edits;
  if (s.alignment.score >= std::numeric_limits<std::int32_t>::min() &&
      s.alignment.score <= std::numeric_limits<std::int32_t>::max()) out << "\tAS:i:" << s.alignment.score;
  if (!s.anchor_ids.empty()) { out << "\tai:Z:"; IdList(out, s.anchor_ids); }
  if (!s.gap_ids.empty()) { out << "\tgi:Z:"; IdList(out, s.gap_ids); }
}
void CheckOutput(std::ostream& out) {
  if (!out) throw AlignmentError("failed to write centromere alignment export");
}
void QueryFragment(std::ostream& out, const SequenceRecord& query, const AlignmentRecord& a) {
  if (a.strand == Strand::Forward)
    out.write(query.bases.data() + a.query_begin, static_cast<std::streamsize>(a.query_end - a.query_begin));
  else for (Position q = a.query_end; q > a.query_begin; ) out.put(Complement(query.bases[--q]));
}
// A map of disjoint query intervals permits O(log n) admission per segment.
bool Admit(std::map<Position, Position>& occupied, Position begin, Position end) {
  const auto next = occupied.lower_bound(begin);
  if (next != occupied.end() && next->first < end) return false;
  if (next != occupied.begin() && std::prev(next)->second > begin) return false;
  occupied.emplace(begin, end);
  return true;
}
void MafRow(std::ostream& out, bool reference_row, const SequenceRecord& reference,
            const SequenceRecord& query, const AlignmentRecord& a) {
  const auto& sequence = reference_row ? reference : query;
  const auto begin = reference_row ? a.reference_begin : a.query_begin;
  const auto end = reference_row ? a.reference_end : a.query_end;
  const bool reverse = !reference_row && a.strand == Strand::Reverse;
  out << "s " << (reference_row ? "ref." : "query.") << sequence.name << ' '
      << (reverse ? sequence.size() - end : begin) << ' ' << end - begin << ' '
      << (reverse ? '-' : '+') << ' ' << sequence.size() << ' ';
  Position r = a.reference_begin, q = 0;
  for (const auto& op : a.cigar) {
    const bool gap = reference_row ? op.operation == 'I' : op.operation == 'D';
    for (Length k = 0; k < op.length; ++k) {
      if (gap) out.put('-');
      else if (reference_row) out.put(reference.bases[r + k]);
      else out.put(reverse ? Complement(query.bases[a.query_end - q - k - 1])
                           : query.bases[a.query_begin + q + k]);
    }
    if (op.operation != 'I') r += op.length;
    if (op.operation != 'D') q += op.length;
  }
  out << '\n';
}
}  // namespace

void WriteCentromereStitchPaf(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  const auto entries = Validate(reference, query, result);
  for (const auto& entry : entries) {
    const auto& a = entry.segment->alignment;
    out << query.name << '\t' << query.size() << '\t' << a.query_begin << '\t' << a.query_end << '\t'
        << Sign(a.strand) << '\t' << reference.name << '\t' << reference.size() << '\t'
        << a.reference_begin << '\t' << a.reference_end << '\t' << entry.counts.matches << '\t'
        << entry.counts.columns << "\t255\tcg:Z:";
    Cigar(out, a.cigar); Tags(out, entry, result); out << '\n';
  }
  CheckOutput(out);
}

void WriteCentromereStitchSam(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  auto entries = Validate(reference, query, result);
  CheckSamNames(reference, query);
  std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
    if (a.segment->path_role != b.segment->path_role) return a.segment->path_role == "primary";
    return a.segment->segment_id < b.segment->segment_id;
  });
  if (!entries.empty()) {
    const auto representative = std::max_element(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
      if (a.segment->path_role != b.segment->path_role) return b.segment->path_role == "primary";
      if (a.counts.paired != b.counts.paired) return a.counts.paired < b.counts.paired;
      return a.segment->segment_id > b.segment->segment_id;
    });
    std::rotate(entries.begin(), representative, std::next(representative));
  }
  out << "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:" << reference.name << "\tLN:" << reference.size()
      << "\n@PG\tID:ramag-centromere-stitch\tPN:RaMA-G\n"
      << "@CO\tSegments only; hard clipping retains full-query coordinates. MAPQ unavailable.\n"
      << "@CO\tpr=primary/copy is algorithm role, not SAM flag. Secondary records retain accepted extra relations; do not filter for full relation evaluation.\n"
      << "@CO\tNo SA tags: path/segment IDs link output without quadratic tag expansion. Status="
      << (result.complete ? "complete" : "partial") << "; complete scope=eligible_internal_gaps.\n";
  if (entries.empty()) {
    out << query.name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << (query.bases.empty() ? "*" : query.bases)
        << "\t*\tst:Z:" << (result.complete ? "complete" : "partial") << '\n';
  }
  std::map<Position, Position> occupied;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    const auto& entry = entries[i];
    const auto& a = entry.segment->alignment;
    const bool canonical = Admit(occupied, a.query_begin, a.query_end);
    const unsigned flag = (a.strand == Strand::Reverse ? 16U : 0U) | (i ? (canonical ? 2048U : 256U) : 0U);
    const Length left = a.strand == Strand::Forward ? a.query_begin : query.size() - a.query_end;
    const Length right = a.strand == Strand::Forward ? query.size() - a.query_end : a.query_begin;
    out << query.name << '\t' << flag << '\t' << reference.name << '\t' << a.reference_begin + 1 << "\t255\t";
    detail::WriteStitchSamCigarRun(out, left, 'H');
    for (const auto& op : a.cigar) detail::WriteStitchSamCigarRun(out, op.length, op.operation);
    detail::WriteStitchSamCigarRun(out, right, 'H');
    out << "\t*\t0\t0\t"; QueryFragment(out, query, a); out << "\t*";
    Tags(out, entry, result); out << '\n';
  }
  CheckOutput(out);
}

void WriteCentromereStitchMaf(std::ostream& out, const SequenceRecord& reference,
                             const SequenceRecord& query, const CentromereStitchResult& result) {
  const auto entries = Validate(reference, query, result);
  out << "##maf version=1 program=RaMA-G scoring=" << CentromereStitchScoreModelName(result.score_model)
      << "\n# ref. and query. prefix original FASTA names to distinguish sequence roles.\n"
      << "# status=" << (result.complete ? "complete" : "partial")
      << " complete_scope=eligible_internal_gaps; independent blocks may share coordinates.\n\n";
  for (const auto& entry : entries) {
    const auto& s = *entry.segment;
    out << "# ramag si=" << s.segment_id << " pi=" << s.path_id << " pr=" << s.path_role
        << " ri=" << s.region_id << " bi=" << s.branch_id
        << " sm=" << CentromereStitchScoreModelName(result.score_model)
        << " st=" << (result.complete ? "complete" : "partial") << " sc=" << s.alignment.score;
    if (!s.anchor_ids.empty()) { out << " ai="; IdList(out, s.anchor_ids); }
    if (!s.gap_ids.empty()) { out << " gi="; IdList(out, s.gap_ids); }
    out << "\na score=" << s.alignment.score << '\n';
    MafRow(out, true, reference, query, s.alignment);
    MafRow(out, false, reference, query, s.alignment);
    out << '\n';
  }
  CheckOutput(out);
}
}  // namespace ramag

namespace ramag {
using namespace centromere_refinement_detail;
namespace {
void WriteEndpointReviews(std::ostream &o,
                          const CentromereEndpointResult &out,bool contexts) {
  o << "review_id\tregion_id\tleft_anchor_id\tanchor_id\tright_anchor_"
       "id\treference_begin\treference_end\tquery_begin\tquery_end\tleft_"
       "gap\tright_gap\told_exact_bases\tnew_exact_bases\tnew_mismatch_"
       "bases\told_score\tnew_"
       "score\tscored\taccepted\treason\tengine\tpeak_kernel_bytes\tdp_"
       "cells\twfa_states\tseconds\told_cigar\tnew_cigar";
  if(contexts)o << "\told_unique_context_starts\tnew_unique_context_starts";
  o << '\n' << std::setprecision(17);
  for (const auto &a : out.reviews) {
    o << a.id << '\t' << a.region_id << '\t' << a.left_anchor_id << '\t'
      << a.anchor_id << '\t' << a.right_anchor_id << '\t' << a.reference_begin
      << '\t' << a.reference_end << '\t' << a.query_begin << '\t' << a.query_end
      << '\t' << a.left_gap << '\t' << a.right_gap << '\t' << a.old_exact_bases
      << '\t' << a.new_exact_bases << '\t' << a.new_mismatch_bases << '\t'
      << a.old_score << '\t' << a.new_score << '\t' << a.scored << '\t'
      << a.accepted << '\t' << a.reason << '\t' << a.engine << '\t'
      << a.peak_kernel_bytes << '\t' << a.dp_cells << '\t' << a.wfa_states
      << '\t' << a.seconds << '\t';
    Cigar(o, a.old_cigar);
    o << '\t';
    Cigar(o, a.new_cigar);
    if(contexts)o << '\t' << a.old_unique_context_starts << '\t' << a.new_unique_context_starts;
    o << '\n';
  }
}
} // namespace
void WriteCentromereEndpointReviews(std::ostream &o,const CentromereEndpointResult &out) {
  WriteEndpointReviews(o,out,false);
}
void WriteCentromereCrossRegionReviews(std::ostream &o,const CentromereEndpointResult &out) {
  WriteEndpointReviews(o,out,true);
}
} // namespace ramag

#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag {
using namespace centromere_refinement_detail;
namespace centromere_refinement_detail {
// Exact 2-bit 32-mers use all 64 bits, so equality has no hash collision
// semantics. Counts saturate at two because only physical uniqueness matters.
int ContextBits(char c) {
  switch (c) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}
ContextCounts CountContexts(const SequenceRecord &sequence,
                            const CentromereEndpointOptions &options) {
  ContextCounts counts;
  std::uint64_t key = 0;
  unsigned run = 0;
  for (std::size_t i = 0; i < sequence.bases.size(); ++i) {
    if ((i % 4096) == 0 && options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-block-context");
    const int bits = ContextBits(sequence.bases[i]);
    if (bits < 0) {
      key = 0;
      run = 0;
      continue;
    }
    key = (key << 2U) | static_cast<std::uint64_t>(bits);
    run = std::min(32U, run + 1);
    if (run == 32U) {
      auto &count = counts[key];
      if (count < 2)
        ++count;
    }
  }
  return counts;
}
bool DistinguishableContext(const A &anchor, const SequenceRecord &reference,
                            const ContextCounts &rc, const ContextCounts &qc) {
  if (anchor.length < 64)
    return false;
  std::uint64_t key = 0;
  Length first_end = 0;
  bool first = false;
  for (Length i = 0; i < anchor.length; ++i) {
    const int bits = ContextBits(reference.bases[anchor.reference_begin + i]);
    // Input anchors have already passed canonical exact-match validation.
    key = (key << 2U) | static_cast<std::uint64_t>(bits);
    if (i < 31)
      continue;
    const auto r = rc.find(key), q = qc.find(key);
    if (r == rc.end() || q == qc.end() || r->second != 1 || q->second != 1)
      continue;
    const Length begin = i - 31;
    if (first && begin >= first_end)
      return true;
    if (!first) {
      first = true;
      first_end = i + 1;
    }
  }
  return false;
}
bool CigarOffset(Position r, Position q, Position rb, Position qb, bool dr,
                 bool dq, Length n, Length &offset) {
  if (rb < r || qb < q)
    return false;
  offset = dr ? rb - r : qb - q;
  return offset <= n && rb - r == (dr ? offset : 0) &&
         qb - q == (dq ? offset : 0);
}
bool SliceStitchedCigar(const CentromereStitchResult &stitched,
                        std::uint64_t region, Position rb, Position re,
                        Position qb, Position qe, std::vector<CigarOp> &out) {
  for (const auto &segment : stitched.segments) {
    const auto &a = segment.alignment;
    if (segment.path_role != "primary" || segment.region_id != region ||
        a.strand != Strand::Forward || a.reference_begin > rb ||
        a.reference_end < re || a.query_begin > qb || a.query_end < qe)
      continue;
    Position r = a.reference_begin, q = a.query_begin;
    bool active = false;
    std::vector<CigarOp> sliced;
    for (auto op : a.cigar) {
      if (op.operation != '=' && op.operation != 'X' && op.operation != 'I' &&
          op.operation != 'D')
        throw AlignmentError("block review requires extended stitched CIGAR");
      const bool dr = op.operation != 'I', dq = op.operation != 'D';
      Length begin = 0, end = 0;
      if (!active) {
        if (!CigarOffset(r, q, rb, qb, dr, dq, op.length, begin)) {
          r += dr ? op.length : 0;
          q += dq ? op.length : 0;
          continue;
        }
        active = true;
      }
      if (CigarOffset(r, q, re, qe, dr, dq, op.length, end) && end >= begin) {
        Emit(sliced, op.operation, end - begin);
        out = std::move(sliced);
        return true;
      }
      Emit(sliced, op.operation, op.length - begin);
      r += dr ? op.length : 0;
      q += dq ? op.length : 0;
    }
  }
  return false;
}
} // namespace centromere_refinement_detail

CentromereEndpointResult
RefineCentromereAnchorBlocks(const SequenceRecord &r, const SequenceRecord &q,
                             const CentromereDirectedAnchorResult &input,
                             const CentromereCopyAnchorResult *copies,
                             const CentromereStitchResult &stitched,
                             const CentromereEndpointOptions &options) {
  if (!input.complete || (copies && !copies->complete) || stitched.aborted ||
      stitched.reference_length != r.size() ||
      stitched.query_length != q.size())
    throw AlignmentError("block review requires valid completed anchor stages");
  if (!options.max_reviews || !std::isfinite(options.stage_seconds_limit) ||
      options.stage_seconds_limit <= 0 ||
      !std::isfinite(options.window_seconds_limit) ||
      options.window_seconds_limit <= 0 ||
      options.kernel.score_model !=
          CentromereStitchScoreModel::Centrolign3Piece ||
      stitched.score_model != options.kernel.score_model)
    throw AlignmentError("invalid block review scoring/resources");
  const auto started = Clock::now();
  const auto check = [&] {
    if (options.kernel.interruption_callback)
      options.kernel.interruption_callback("centromere-block-review");
  };
  check();
  CentromereEndpointResult out;
  out.anchors = input;
  std::map<std::uint64_t, std::vector<const A *>> paths;
  std::set<std::uint64_t> ids;
  for (const auto &a : input.final) {
    CheckAnchor(a, r, q);
    if (!ids.insert(a.id).second)
      throw AlignmentError("duplicate block review anchor ID");
    paths[a.region_id].push_back(&a);
  }
  for (auto &[id, path] : paths) {
    (void)id;
    std::sort(path.begin(), path.end(),
              [](auto a, auto b) { return a->chain_index < b->chain_index; });
    for (std::size_t i = 0; i < path.size(); ++i) {
      const auto *a = path[i];
      if (a->chain_index != i + 1)
        throw AlignmentError("block review chain index");
      if (i) {
        const auto *b = path[i - 1];
        if (a->strand != b->strand ||
            b->reference_begin + b->length > a->reference_begin ||
            (a->strand == Strand::Forward
                 ? b->query_begin + b->length > a->query_begin
                 : a->query_begin + a->length > b->query_begin))
          throw AlignmentError("block review chain order");
      }
    }
  }
  const auto rc = CountContexts(r, options), qc = CountContexts(q, options);
  std::vector<Replacement> accepted;
  for (const auto &[region, path] : paths) {
    if (path.empty() || path.front()->strand != Strand::Forward)
      continue;
    std::vector<std::size_t> strong;
    for (std::size_t i = 0; i < path.size(); ++i) {
      check();
      if (DistinguishableContext(*path[i], r, rc, qc))
        strong.push_back(i);
    }
    for (std::size_t k = 1; k < strong.size(); ++k) {
      check();
      const std::size_t li = strong[k - 1], zi = strong[k];
      // A single-root window is already handled by the original stage.
      if (zi - li < 3)
        continue;
      const auto &l = *path[li], &z = *path[zi];
      CentromereEndpointReview review;
      review.id = out.reviews.size() + 1;
      review.region_id = region;
      review.left_anchor_id = l.id;
      review.anchor_id = path[li + 1]->id;
      review.right_anchor_id = z.id;
      review.reference_begin = l.reference_begin + l.length;
      review.reference_end = z.reference_begin;
      review.query_begin = l.query_begin + l.length;
      review.query_end = z.query_begin;
      review.left_gap =
          std::max(path[li + 1]->reference_begin - review.reference_begin,
                   path[li + 1]->query_begin - review.query_begin);
      review.right_gap = std::max(
          review.reference_end - path[zi - 1]->reference_begin -
              path[zi - 1]->length,
          review.query_end - path[zi - 1]->query_begin - path[zi - 1]->length);
      const Length rn = review.reference_end - review.reference_begin,
                   qn = review.query_end - review.query_begin;
      if (zi - li - 1 > 256 || rn > 65536 || qn > 65536)
        review.reason = "block_geometry_budget";
      bool has_root = false;
      for (const auto &a : input.initial)
        if (a.region_id == region && a.strand == Strand::Forward &&
            a.reference_begin >= review.reference_begin &&
            a.reference_begin + a.length <= review.reference_end &&
            a.query_begin >= review.query_begin &&
            a.query_begin + a.length <= review.query_end) {
          has_root = true;
          break;
        }
      if (review.reason.empty() && !has_root)
        review.reason = "no_root_in_block";
      if (review.reason.empty() && copies)
        for (const auto &a : copies->anchors)
          if (Touch(a.reference_begin, a.reference_begin + a.length,
                    review.reference_begin, review.reference_end) ||
              Touch(a.query_begin, a.query_begin + a.length, review.query_begin,
                    review.query_end)) {
            review.reason = "accepted_copy_footprint";
            break;
          }
      if (review.reason.empty())
        for (const auto &a : accepted)
          if (Touch(a.rb, a.re, review.reference_begin, review.reference_end) ||
              Touch(a.qb, a.qe, review.query_begin, review.query_end)) {
            review.reason = "overlapping_accepted_window";
            break;
          }
      if (review.reason.empty() &&
          !SliceStitchedCigar(stitched, region, review.reference_begin,
                              review.reference_end, review.query_begin,
                              review.query_end, review.old_cigar))
        review.reason = "unresolved_or_missing_old_path";
      const auto rv = std::string_view(r.bases).substr(review.reference_begin,
                                                       rn),
                 qv = std::string_view(q.bases).substr(review.query_begin, qn);
      if (review.reason.empty()) {
        review.old_score = detail::ScoreCentromereStitchCigar(
            rv, qv, review.old_cigar, options.kernel.score_model);
        for (auto op : review.old_cigar)
          if (op.operation == '=')
            review.old_exact_bases += op.length;
        if (review.old_exact_bases >= std::min(rn, qn))
          review.reason = "no_exact_headroom";
        else if (rv == qv)
          review.reason = "identical_window_ancestry_ambiguous";
      }
      if (review.reason.empty() &&
          (out.kernel_calls >= options.max_reviews ||
           Seconds(started) >= options.stage_seconds_limit))
        review.reason = "review_budget";
      if (!review.reason.empty()) {
        out.reviews.push_back(std::move(review));
        continue;
      }
      auto kernel = options.kernel;
      kernel.progress_callback = {};
      kernel.observation_callback = {};
      kernel.gap_seconds_limit =
          std::min({kernel.gap_seconds_limit, options.window_seconds_limit,
                    options.stage_seconds_limit - Seconds(started)});
      if (kernel.gap_seconds_limit <= 0) {
        review.reason = "review_budget";
        out.reviews.push_back(std::move(review));
        continue;
      }
      const auto begin = Clock::now();
      ++out.kernel_calls;
      try {
        auto proposal = detail::AlignCentromereBlockDp(rv, qv, kernel);
        review.new_score = detail::ScoreCentromereStitchCigar(
            rv, qv, proposal.cigar, kernel.score_model);
        if (review.new_score != proposal.score)
          throw AlignmentError("block kernel score/replay differs");
        review.scored = true;
        review.engine = proposal.engine;
        review.peak_kernel_bytes = proposal.peak_kernel_bytes;
        review.dp_cells = proposal.dp_cells;
        review.wfa_states = proposal.wfa_states;
        review.new_cigar = std::move(proposal.cigar);
        Replacement p{region,
                      review.reference_begin,
                      review.reference_end,
                      review.query_begin,
                      review.query_end,
                      {}};
        Position rb = p.rb, qb = p.qb;
        for (auto op : review.new_cigar) {
          if (op.operation == 'X')
            review.new_mismatch_bases += op.length;
          if (op.operation == '=') {
            p.exact.push_back({0, region, 0, rb, qb, op.length, Strand::Forward,
                               static_cast<double>(op.length)});
            review.new_exact_bases += op.length;
          }
          if (op.operation != 'I')
            rb += op.length;
          if (op.operation != 'D')
            qb += op.length;
        }
        if (review.new_mismatch_bases)
          review.reason = "new_substitutions_preserved";
        else if (review.new_score <= review.old_score)
          review.reason = review.new_score == review.old_score
                              ? "equal_score_preserved"
                              : "lower_score_preserved";
        else if (review.new_exact_bases <= review.old_exact_bases)
          review.reason = "no_exact_coverage_gain";
        else {
          review.accepted = true;
          review.reason = "strict_score_and_exact_gain";
          accepted.push_back(std::move(p));
          ++out.accepted;
        }
      } catch (const detail::StitchGapResourceLimit &error) {
        review.reason = error.reason;
        review.engine = error.engine;
        review.peak_kernel_bytes = error.peak_kernel_bytes;
        review.dp_cells = error.dp_cells;
        review.wfa_states = error.wfa_states;
      }
      review.seconds = Seconds(begin);
      out.reviews.push_back(std::move(review));
    }
  }
  if (!accepted.empty()) {
    Replace(out.anchors.initial, accepted, input.regions);
    Replace(out.anchors.final, accepted, input.regions);
    Refresh(out.anchors);
  }
  check();
  out.seconds = Seconds(started);
  return out;
}
} // namespace ramag

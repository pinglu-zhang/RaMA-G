#include "ramag/centromere_direction.hpp"
#include <iomanip>
#include <ostream>

namespace ramag {
namespace {
char Sign(Strand s) {return s==Strand::Forward ? '+' : '-';}
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
     <<Sign(a.strand)<<'\t'<<a.length<<'\t'<<a.source_weight<<'\t'<<a.full_reference_count<<'\t'<<a.full_query_count<<'\t'
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
    o<<'\t'<<Sign(a.strand)<<'\t'<<Sign(a.background)<<'\t'<<a.classification<<'\t'<<a.reason<<'\t'
     <<a.support<<'\t'<<a.opposite_support<<'\t'<<a.context_background_support<<'\t'<<a.competing<<'\t'<<a.left_flank<<'\t'<<a.right_flank;
    auto bounds=[&](bool known,const auto& interval){if(known)o<<'\t'<<interval[0]<<'\t'<<interval[1];else o<<"\t.\t.";};
    bounds(a.left_bounded,a.reference_left);bounds(a.right_bounded,a.reference_right);
    bounds(a.left_bounded,a.query_left);bounds(a.right_bounded,a.query_right);
    o<<'\t';List(o,a.anchors);o<<'\t'<<a.independent_segments<<'\t'<<a.longest_match<<'\t';List(o,a.enclosing_opposite_candidates);o<<'\n';
  }
}
} // namespace ramag

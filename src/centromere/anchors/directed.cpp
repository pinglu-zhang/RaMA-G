#include "ramag/centromere_directed_anchors.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <numeric>
#include <set>
#include <string_view>
#include <utility>

namespace ramag {
namespace {
using Clock = std::chrono::steady_clock;
using Range = std::pair<Position, Position>;
char Sign(Strand s) { return s == Strand::Forward ? '+' : '-'; }
bool Canonical(char c) { return c=='A' || c=='C' || c=='G' || c=='T'; }
char Complement(char c) {
  switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C';
    case 'T': return 'A'; default: return 'N'; }
}
std::string ReverseComplement(std::string s) {
  std::reverse(s.begin(),s.end());
  for (char& c:s) c=Complement(c);
  return s;
}
void Check(const CentromereAlignmentOptions& o) {
  if (o.interruption_callback) o.interruption_callback("centromere-directed-anchors");
}
std::vector<std::string_view> Fields(std::string_view line) {
  std::vector<std::string_view> out;
  while (true) {
    const auto p=line.find('\t');
    out.push_back(line.substr(0,p));
    if (p==std::string_view::npos) return out;
    line.remove_prefix(p+1);
  }
}
std::uint64_t Number(std::string_view s,std::size_t line) {
  std::uint64_t n{};
  const auto [end,error]=std::from_chars(s.data(),s.data()+s.size(),n);
  if (error!=std::errc{} || end!=s.data()+s.size() || s.empty())
    throw CentromereInputError("invalid unsigned coordinate/id at direction line "+std::to_string(line));
  return n;
}
void ValidateRegions(std::span<const CentromereDirectedRegion> regions,Length query_length) {
  std::set<std::uint64_t> ids;
  std::vector<Range> bounds;
  for (const auto& r:regions) {
    if (!r.id || !ids.insert(r.id).second)
      throw CentromereInputError("direction region IDs must be positive and unique");
    if (r.query_begin>=r.query_end || r.query_end>query_length)
      throw CentromereInputError("direction region is empty or outside query");
    if (r.strand!=Strand::Forward && r.strand!=Strand::Reverse)
      throw CentromereInputError("invalid direction region strand");
    bounds.emplace_back(r.query_begin,r.query_end);
  }
  std::sort(bounds.begin(),bounds.end());
  Position end=0;
  for (const auto& r:bounds) {
    if (r.first<end) throw CentromereInputError("direction regions overlap or repeat");
    end=r.second;
  }
}
std::vector<Range> Union(std::vector<Range> ranges) {
  std::sort(ranges.begin(),ranges.end());
  std::vector<Range> result;
  for (const auto& r:ranges) {
    if (r.first==r.second) continue;
    if (result.empty() || result.back().second<r.first) result.push_back(r);
    else result.back().second=std::max(result.back().second,r.second);
  }
  return result;
}
Length Bases(const std::vector<Range>& ranges) {
  Length n=0;
  for (const auto& r:ranges) n+=r.second-r.first;
  return n;
}
std::vector<Range> ComplementRanges(Position begin,Position end,const std::vector<Range>& ranges) {
  std::vector<Range> out;
  for (const auto& r:ranges) {
    if (begin<r.first) out.emplace_back(begin,r.first);
    begin=r.second;
  }
  if (begin<end) out.emplace_back(begin,end);
  return out;
}
Position GlobalQuery(const CentromereDirectedRegion& r,Position local,Length length) {
  if (local>r.query_end-r.query_begin || length>r.query_end-r.query_begin-local)
    throw AlignmentError("directed anchor/trace outside query segment");
  return r.strand==Strand::Forward ? r.query_begin+local : r.query_end-local-length;
}
void AppendAnchors(std::span<const CentromereAnchor> input,
    const CentromereDirectedRegion& region,const SequenceRecord& reference,
    const SequenceRecord& query,std::vector<CentromereDirectedAnchor>& output) {
  Position previous_reference=0,previous_query=0;
  std::uint64_t chain_index=0;
  for (const auto& a:input) {
    if (!a.length || a.reference_begin>reference.size() ||
        a.length>reference.size()-a.reference_begin ||
        a.reference_begin<previous_reference || a.query_begin<previous_query ||
        !std::isfinite(a.weight)) throw AlignmentError("invalid directed anchor order, bounds or weight");
    const auto qb=GlobalQuery(region,a.query_begin,a.length);
    for (Length i=0;i<a.length;++i) {
      const char r=reference.bases[a.reference_begin+i];
      const char q=region.strand==Strand::Forward ? query.bases[qb+i]
        : Complement(query.bases[qb+a.length-1-i]);
      if (!Canonical(r) || r!=q) throw AlignmentError("directed anchor is not canonical exact match");
    }
    output.push_back({static_cast<std::uint64_t>(output.size()+1),region.id,++chain_index,
                      a.reference_begin,qb,a.length,region.strand,a.weight});
    previous_reference=a.reference_begin+a.length;
    previous_query=a.query_begin+a.length;
  }
}
void RequireInitialRetained(std::span<const CentromereAnchor> initial,
                            std::span<const CentromereAnchor> final) {
  std::size_t j=0;
  for (const auto& a:initial) {
    while (j<final.size() && final[j].reference_begin+final[j].length<=a.reference_begin) ++j;
    if (j==final.size()) throw AlignmentError("recursive result lost an initial anchor");
    const auto& b=final[j];
    if (a.reference_begin<b.reference_begin || a.query_begin<b.query_begin ||
        a.reference_begin-b.reference_begin!=a.query_begin-b.query_begin ||
        a.length>b.length || a.reference_begin-b.reference_begin>b.length-a.length)
      throw AlignmentError("recursive result changed an initial anchor correspondence");
  }
}
std::vector<Range> AnchorRanges(std::span<const CentromereAnchor> anchors,bool reference) {
  std::vector<Range> out;
  out.reserve(anchors.size());
  for (const auto& a:anchors) {
    const auto b=reference ? a.reference_begin : a.query_begin;
    out.emplace_back(b,b+a.length);
  }
  return Union(std::move(out));
}
void WriteAnchors(std::ostream& out,std::span<const CentromereDirectedAnchor> anchors) {
  out<<"anchor_id\tregion_id\tchain_index\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tlength\tweight\n";
  out<<std::setprecision(17);
  for (const auto& a:anchors)
    out<<a.id<<'\t'<<a.region_id<<'\t'<<a.chain_index<<'\t'<<a.reference_begin<<'\t'
       <<a.reference_begin+a.length<<'\t'<<a.query_begin<<'\t'<<a.query_begin+a.length<<'\t'
       <<Sign(a.strand)<<'\t'<<a.length<<'\t'<<a.weight<<'\n';
}
} // namespace

std::vector<CentromereDirectedRegion> ReadCentromereDirectedRegions(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw CentromereInputError("cannot open fixed direction file: "+path.string());
  std::string line;
  std::size_t number=0;
  bool header=false,canonical=false;
  std::vector<CentromereDirectedRegion> result;
  while (std::getline(input,line)) {
    ++number;
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line.empty() || line.front()=='#') continue;
    if (!header) {
      if (line=="block_id\tquery_begin\tquery_end\tstrand") canonical=false;
      else if (line=="region_id\tsource_line\tquery_begin\tquery_end\tstrand") canonical=true;
      else throw CentromereInputError("unsupported fixed direction TSV header");
      header=true;
      continue;
    }
    const auto f=Fields(line);
    if (f.size()!=(canonical ? 5U : 4U))
      throw CentromereInputError("wrong fixed direction column count at line "+std::to_string(number));
    const std::size_t offset=canonical ? 2 : 1;
    (void)Number(f[0],number); // Legacy block ID is deliberately not a pairing key.
    const auto strand=f[offset+2];
    if (strand!="+" && strand!="-")
      throw CentromereInputError("invalid fixed direction strand at line "+std::to_string(number));
    const auto source=canonical ? Number(f[1],number) : static_cast<std::uint64_t>(number);
    if (source>std::numeric_limits<std::size_t>::max()) throw CentromereInputError("direction source line overflows");
    result.push_back({canonical ? Number(f[0],number) : static_cast<std::uint64_t>(result.size()+1),
      Number(f[offset],number),Number(f[offset+1],number),
      strand=="+" ? Strand::Forward : Strand::Reverse,static_cast<std::size_t>(source)});
  }
  if (input.bad()) throw CentromereInputError("cannot read fixed direction file: "+path.string());
  if (!header) throw CentromereInputError("fixed direction TSV is missing its header");
  ValidateRegions(result,std::numeric_limits<Length>::max());
  return result;
}

CentromereDirectedAnchorResult FindCentromereAnchorsInDirections(
    const SequenceRecord& reference,const SequenceRecord& query,
    std::span<const CentromereDirectedRegion> regions,const CentromereAlignmentOptions& options) {
  if (!CentromereBackendAvailable())
    throw DependencyError("centromere backend not built; enable RAMAG_ENABLE_CENTROMERE on Linux");
  if (!options.worker_threads) throw CentromereInputError("centromere threads must be positive");
  if (reference.bases.empty() || query.bases.empty()) throw CentromereInputError("directed anchor sequences must be nonempty");
  for (const auto* s:{&reference,&query})
    for (const char c:s->bases)
      if (!Canonical(c) && c!='N') throw CentromereInputError("directed anchors require normalized A/C/G/T/N");
  ValidateRegions(regions,query.size());
  Check(options);
  const auto start=Clock::now();
  CentromereDirectedAnchorResult result;
  result.regions.assign(regions.begin(),regions.end());
  std::vector<Range> known;
  for (const auto& region:regions) known.emplace_back(region.query_begin,region.query_end);
  for (const auto& [b,e]:ComplementRanges(0,query.size(),Union(std::move(known))))
    result.unknown.push_back({0,b,e,Strand::Forward});
  std::vector<std::vector<Range>> reference_ranges;
  reference_ranges.reserve(regions.size());
  // Preserve each unconstrained search as evidence. A flanking window is a
  // repair only for a root chain that contradicts that window, not a reason to
  // change rarity and seed a different recursion in an already localized chain.
  struct CachedSearch {
    CentromereAnchorResult anchors;
    std::vector<CentromereRecursionNode> trace;
    double seconds{};
  };
  std::vector<std::optional<CachedSearch>> baseline(regions.size());
  if (options.progress_callback) options.progress_callback("centromere-directed-regions",0,regions.size());
  for (std::size_t i=0;i<regions.size();++i) {
    const auto& region=regions[i];
    Check(options);
    const auto begin=Clock::now();
    auto& cached=baseline[i].emplace();
    SequenceRecord segment{query.numeric_id,query.name,query.header,
      query.bases.substr(region.query_begin,region.query_end-region.query_begin)};
    if (region.strand==Strand::Reverse) segment.bases=ReverseComplement(std::move(segment.bases));
    auto local_options=options;
    local_options.corrected_recursive_boundaries=true;
    local_options.recursion_observer=[&](const CentromereRecursionNode& node) {
      cached.trace.push_back(node);
    };
    cached.anchors=FindCentromereAnchors(reference,segment,local_options);
    cached.seconds=std::chrono::duration<double>(Clock::now()-begin).count();
  }
  std::vector<std::size_t> order(regions.size()),rank(regions.size());
  std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[&](auto a,auto b) {
    return regions[a].query_begin<regions[b].query_begin;
  });
  for (std::size_t i=0;i<order.size();++i) rank[order[i]]=i;
  for (std::size_t region_index=0;region_index<regions.size();++region_index) {
    const auto& region=regions[region_index];
    Check(options);
    const auto task_start=Clock::now();
    SequenceRecord oriented{query.numeric_id,query.name,query.header,
      query.bases.substr(region.query_begin,region.query_end-region.query_begin)};
    if (region.strand==Strand::Reverse) oriented.bases=ReverseComplement(std::move(oriented.bases));
    CentromereDirectedTask task;
    task.region=region;task.reference_length=reference.size();
    task.coverage_guard=baseline[region_index]->anchors.coverage_guard;
    task.reference_end=reference.size();
    // A reverse segment bracketed by two forward segments is a local inversion
    // candidate. Ordered flanking correspondences delimit its reference search;
    // otherwise retain the unconstrained route (including duplicated targets).
    const auto position=rank[region_index];
    if (region.strand==Strand::Reverse && position>0 && position+1<order.size()) {
      const auto left=order[position-1],right=order[position+1];
      if (regions[left].strand==Strand::Forward && regions[right].strand==Strand::Forward &&
          !baseline[left]->anchors.final.empty() && !baseline[right]->anchors.final.empty()) {
        const auto& a=baseline[left]->anchors.final.back();
        const auto& b=baseline[right]->anchors.final.front();
        const auto begin=a.reference_begin+a.length,end=b.reference_begin;
        // The flank proposal must be able to contain the complete query span.
        // A narrower box is contradictory evidence (or an indel-rich locus),
        // not permission to discard existing correspondences with a hard cut.
        const auto& root=baseline[region_index]->anchors.initial;
        const bool escapes=std::any_of(root.begin(),root.end(),[&](const auto& anchor) {
          return anchor.reference_begin<begin || anchor.reference_begin+anchor.length>end;
        });
        if (begin<end && end-begin>=region.query_end-region.query_begin && escapes) {
          task.reference_begin=begin;task.reference_end=end;
        }
      }
    }
    auto local_options=options;
    local_options.corrected_recursive_boundaries=true;
    local_options.recursion_observer=[&](const CentromereRecursionNode& node) {
      if (node.reference_begin>node.reference_end || node.reference_end>reference.size() ||
          node.query_begin>node.query_end || node.query_end>oriented.size())
        throw AlignmentError("directed recursion trace outside input");
      result.trace.push_back({region.id,node});
      ++task.recursion_node_count;
      task.max_recursion_depth=std::max(task.max_recursion_depth,node.depth);
      if (options.recursion_observer) options.recursion_observer(node);
    };
    CentromereAnchorResult anchors;
    AppendAnchors(baseline[region_index]->anchors.initial,region,reference,query,result.full_reference_initial);
    if (task.reference_begin==0 && task.reference_end==reference.size()) {
      for (const auto& node:baseline[region_index]->trace) local_options.recursion_observer(node);
      anchors=baseline[region_index]->anchors;
    } else {
      SequenceRecord target{reference.numeric_id,reference.name,reference.header,
        reference.bases.substr(task.reference_begin,task.reference_end-task.reference_begin)};
      const auto observe=local_options.recursion_observer;
      local_options.recursion_observer=[&](CentromereRecursionNode node) {
        node.reference_begin+=task.reference_begin;node.reference_end+=task.reference_begin;
        observe(node);
      };
      anchors=FindCentromereAnchors(target,oriented,local_options);
      task.window_coverage_guard=anchors.coverage_guard;
      for (auto* list:{&anchors.initial,&anchors.final})
        for (auto& anchor:*list) anchor.reference_begin+=task.reference_begin;
    }
    AppendAnchors(anchors.initial,region,reference,query,result.initial);
    AppendAnchors(anchors.final,region,reference,query,result.final);
    RequireInitialRetained(anchors.initial,anchors.final);
    task.initial_anchor_count=anchors.initial.size();task.final_anchor_count=anchors.final.size();
    task.initial_reference_bases=Bases(AnchorRanges(anchors.initial,true));
    task.initial_query_bases=Bases(AnchorRanges(anchors.initial,false));
    const auto final_ref=AnchorRanges(anchors.final,true);
    const auto final_query=AnchorRanges(anchors.final,false);
    task.final_reference_bases=Bases(final_ref);task.final_query_bases=Bases(final_query);
    if (task.final_query_bases<task.initial_query_bases) throw AlignmentError("recursive coverage decreased");
    task.recursive_new_query_bases=task.final_query_bases-task.initial_query_bases;
    reference_ranges.push_back(final_ref);
    for (const auto& [b,e]:ComplementRanges(0,oriented.size(),final_query)) {
      const auto qb=GlobalQuery(region,b,e-b);
      result.uncovered.push_back({region.id,qb,qb+e-b,region.strand});
      task.unanchored_query_bases+=e-b;
      ++task.unanchored_interval_count;
      task.longest_unanchored_query_bases=std::max(task.longest_unanchored_query_bases,e-b);
    }
    task.status=anchors.final.empty() ? "no_anchors" : "complete";
    task.total_seconds=std::chrono::duration<double>(Clock::now()-task_start).count();
    task.total_seconds+=baseline[region_index]->seconds;
    result.tasks.push_back(std::move(task));
    Check(options);
    if (options.progress_callback) options.progress_callback("centromere-directed-regions",result.tasks.size(),regions.size());
  }
  // Independent segments may select the same reference copy. Report overlap;
  // do not impose a new cross-segment one-to-one assignment or discard anchors.
  for (std::size_t a=0;a<regions.size();++a) {
    Check(options);
    for (std::size_t b=a+1;b<regions.size();++b) {
      std::size_t i=0,j=0;
      const auto& x=reference_ranges[a];const auto& y=reference_ranges[b];
      while (i<x.size() && j<y.size()) {
        const auto begin=std::max(x[i].first,y[j].first),end=std::min(x[i].second,y[j].second);
        if (begin<end) result.reference_overlaps.push_back({regions[a].id,regions[b].id,begin,end});
        if (x[i].second<y[j].second) ++i; else ++j;
      }
    }
  }
  std::sort(result.uncovered.begin(),result.uncovered.end(),[](const auto& a,const auto& b) {
    return a.query_begin<b.query_begin;
  });
  Check(options);
  result.total_seconds=std::chrono::duration<double>(Clock::now()-start).count();
  result.complete=true;
  return result;
}

void WriteCentromereDirectedRegions(std::ostream& out,const CentromereDirectedAnchorResult& r) {
  out<<"region_id\tsource_line\tquery_begin\tquery_end\tstrand\n";
  for (const auto& x:r.regions) out<<x.id<<'\t'<<x.source_line<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<Sign(x.strand)<<'\n';
}
void WriteCentromereDirectedInitialAnchors(std::ostream& out,const CentromereDirectedAnchorResult& r) { WriteAnchors(out,r.initial); }
void WriteCentromereDirectedFullReferenceInitialAnchors(std::ostream& out,const CentromereDirectedAnchorResult& r) { WriteAnchors(out,r.full_reference_initial); }
void WriteCentromereDirectedFinalAnchors(std::ostream& out,const CentromereDirectedAnchorResult& r) { WriteAnchors(out,r.final); }
void WriteCentromereDirectedTrace(std::ostream& out,const CentromereDirectedAnchorResult& r) {
  const bool guarded=std::any_of(r.trace.begin(),r.trace.end(),[](const auto& x){return x.node.guard_status!="disabled";});
  out<<"region_id\tnode_id\tparent_id\tdepth\treference_begin\treference_end\toriented_query_begin\toriented_query_end\tquery_begin\tquery_end\tstrand\trare_level\tcandidate_count\tselected_anchor_count\tstop_reason";
  if(guarded)out<<"\tguard_status\tquery_added\treference_added\tcache_bytes";
  out<<'\n';
  std::map<std::uint64_t,CentromereDirectedRegion> regions;
  for (const auto& region:r.regions) regions.emplace(region.id,region);
  for (const auto& x:r.trace) {
    const auto& n=x.node;const auto& region=regions.at(x.region_id);
    const auto qb=GlobalQuery(region,n.query_begin,n.query_end-n.query_begin);
    out<<x.region_id<<'\t'<<n.node_id<<'\t'<<n.parent_id<<'\t'<<n.depth<<'\t'<<n.reference_begin<<'\t'<<n.reference_end
       <<'\t'<<n.query_begin<<'\t'<<n.query_end<<'\t'<<qb<<'\t'<<qb+n.query_end-n.query_begin<<'\t'<<Sign(region.strand)
       <<'\t'<<n.rare_level<<'\t'<<n.candidate_count<<'\t'<<n.selected_anchor_count<<'\t'<<n.stop_reason;
    if(guarded)out<<'\t'<<n.guard_status<<'\t'<<n.query_added<<'\t'<<n.reference_added<<'\t'<<n.cache_bytes;
    out<<'\n';
  }
}
void WriteCentromereDirectedUncovered(std::ostream& out,const CentromereDirectedAnchorResult& r) {
  out<<"region_id\tquery_begin\tquery_end\tstrand\tstatus\n";
  std::vector<std::pair<CentromereDirectedQueryInterval,bool>> intervals;
  for (const auto& x:r.uncovered) intervals.emplace_back(x,false);
  for (const auto& x:r.unknown) intervals.emplace_back(x,true);
  std::sort(intervals.begin(),intervals.end(),[](const auto& a,const auto& b) {return a.first.query_begin<b.first.query_begin;});
  for (const auto& [x,unknown]:intervals)
    out<<x.region_id<<'\t'<<x.query_begin<<'\t'<<x.query_end<<'\t'<<(unknown ? '.' : Sign(x.strand))
       <<'\t'<<(unknown ? "unknown" : "known_unanchored")<<'\n';
}
void WriteCentromereDirectedReferenceOverlaps(std::ostream& out,const CentromereDirectedAnchorResult& r) {
  out<<"first_region_id\tsecond_region_id\treference_begin\treference_end\tlength\n";
  for (const auto& x:r.reference_overlaps)
    out<<x.first_region_id<<'\t'<<x.second_region_id<<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.reference_end-x.reference_begin<<'\n';
}
void WriteCentromereDirectedSummary(std::ostream& out,const CentromereDirectedAnchorResult& r) {
  const bool guarded=std::any_of(r.tasks.begin(),r.tasks.end(),[](const auto& x){return x.coverage_guard.enabled||x.window_coverage_guard.enabled;});
  out<<"region_id\treference_begin\treference_end\tquery_begin\tquery_end\tstrand\tinitial_anchor_count\tfinal_anchor_count\tinitial_reference_bases\tinitial_query_bases\tfinal_reference_bases\tfinal_query_bases\trecursive_new_query_bases\tunanchored_query_bases\tunanchored_interval_count\tlongest_unanchored_query_bases\trecursion_node_count\tmax_recursion_depth\ttotal_seconds\tstatus";
  if(guarded)for(const auto prefix:{"guard_full_","guard_window_"})
    for(const auto name:{"enabled","proposal_seconds_sum","selection_seconds","refill_seconds","baseline_nodes","cached_nodes","cache_exhausted_nodes","proposal_budget_exhausted_nodes","accepted_nodes","reused_subtrees","refill_searches","refill_nodes","refill_nonempty_nodes","cache_bytes_retained"})out<<'\t'<<prefix<<name;
  out<<'\n'<<std::setprecision(17);
  const auto extra=[&](const CentromereCoverageGuardStatistics& g) {
    out<<'\t'<<g.enabled<<'\t'<<g.proposal_seconds_sum<<'\t'<<g.selection_seconds<<'\t'<<g.refill_seconds
       <<'\t'<<g.baseline_nodes<<'\t'<<g.cached_nodes<<'\t'<<g.cache_exhausted_nodes<<'\t'<<g.proposal_budget_exhausted_nodes
       <<'\t'<<g.accepted_nodes<<'\t'<<g.reused_subtrees<<'\t'<<g.refill_searches<<'\t'<<g.refill_nodes
       <<'\t'<<g.refill_nonempty_nodes<<'\t'<<g.cache_bytes_retained;
  };
  for(const auto& x:r.tasks) {
    out<<x.region.id<<'\t'<<x.reference_begin<<'\t'<<x.reference_end<<'\t'<<x.region.query_begin<<'\t'<<x.region.query_end<<'\t'<<Sign(x.region.strand)
       <<'\t'<<x.initial_anchor_count<<'\t'<<x.final_anchor_count<<'\t'<<x.initial_reference_bases<<'\t'<<x.initial_query_bases
       <<'\t'<<x.final_reference_bases<<'\t'<<x.final_query_bases<<'\t'<<x.recursive_new_query_bases<<'\t'<<x.unanchored_query_bases
       <<'\t'<<x.unanchored_interval_count<<'\t'<<x.longest_unanchored_query_bases<<'\t'<<x.recursion_node_count
       <<'\t'<<x.max_recursion_depth<<'\t'<<x.total_seconds<<'\t'<<x.status;
    if(guarded){extra(x.coverage_guard);extra(x.window_coverage_guard);}
    out<<'\n';
  }
}
} // namespace ramag

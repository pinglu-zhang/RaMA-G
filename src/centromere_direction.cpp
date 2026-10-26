#include "ramag/centromere_direction.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <sstream>
#include "centromere_direction_storage.hpp"
#if RAMAG_HAVE_CENTROMERE
#include "anchor.h"
#endif

namespace ramag {
bool CentromereCliDirections() noexcept { return RAMAG_CENTROMERE_DIRECTIONS != 0; }
DirectionSearchStrategy CentromereCliDirectionStrategy() noexcept {
#if RAMAG_DIRECTION_LEGACY_SEARCH
  return DirectionSearchStrategy::LegacyMultiscale;
#else
  return DirectionSearchStrategy::RamaFirstHitFocused;
#endif
}
std::string_view DirectionStrategyName(DirectionSearchStrategy strategy) noexcept {
  return strategy==DirectionSearchStrategy::RamaFirstHitFocused ?
    "rama-first-hit-focused-v1" : "legacy-multiscale-v1";
}
#if RAMAG_HAVE_CENTROMERE
namespace {
using Clock=std::chrono::steady_clock;
void Check(const CentromereDirectionOptions& o) {
  if(o.interruption_callback) o.interruption_callback("centromere-directions");
}
struct MemoryGuard {
  const CentromereDirectionOptions& options;
  std::uint64_t retained{};
  std::uint64_t peak{};
  void Add(std::uint64_t n) {
    if(n>options.evidence_memory_limit_bytes || retained>options.evidence_memory_limit_bytes-n)
      throw AlignmentError("direction evidence memory guard exceeded; retained_estimate_bytes="+
        std::to_string(retained)+" requested_bytes="+std::to_string(n)+" limit_bytes="+
        std::to_string(options.evidence_memory_limit_bytes)+"; analysis incomplete");
    retained+=n;
    peak=std::max(peak,retained);
  }
  void Release(std::uint64_t n) {
    if(n>retained)throw AlignmentError("direction memory accounting underflow");
    retained-=n;
  }
};
Position U(const DirectionAnchor& a, Length q) {
  return a.strand==Strand::Forward ? a.query_begin : q-a.query_begin-a.length;
}
DirectionBox NodeBox(const DirectionChainNode& n, Strand strand, Length q) {
  const auto qb=strand==Strand::Forward ? n.oriented_query_begin : q-n.oriented_query_begin-n.length;
  return {n.reference_begin,n.reference_begin+n.length,qb,qb+n.length};
}
bool Inside(const DirectionBox& a,const DirectionBox& b) {
  return a.reference_begin>=b.reference_begin && a.reference_end<=b.reference_end &&
         a.query_begin>=b.query_begin && a.query_end<=b.query_end;
}
bool Nonempty(const DirectionBox& b) {
  return b.reference_begin<b.reference_end && b.query_begin<b.query_end;
}
std::vector<std::size_t> Trace(const DirectionChain& c,std::size_t endpoint) {
  std::vector<std::size_t> path;
  for(auto i=endpoint;i!=kNoDirectionNode;i=c.nodes[i].predecessor) path.push_back(i);
  std::reverse(path.begin(),path.end()); return path;
}
// Enumerate cells touched by the actual paired diagonal, not the Cartesian
// product of its reference and query spans. Windows use half-width strides.
std::vector<direction_detail::CompactBox> Windows(std::span<const direction_detail::CoveredInterval> mapping,
    const DirectionBox& parent, Length qlength, MemoryGuard& guard) {
  direction_detail::WindowSet boxes;boxes.Insert(parent,guard);
  const auto extent=std::max(parent.reference_end-parent.reference_begin,
                             parent.query_end-parent.query_begin);
  for(Length width=80;width<extent;) {
    const auto step=width/2;
    for(const auto& interval:mapping) {
      const auto a=direction_detail::Mapping(interval,qlength);
      const auto ob=a.strand==Strand::Forward ? parent.query_begin : qlength-parent.query_end;
      const auto oe=a.strand==Strand::Forward ? parent.query_end : qlength-parent.query_begin;
      const auto u=U(a,qlength);
      if(a.reference_begin>=parent.reference_end || u>=oe) continue;
      const auto start=std::max(parent.reference_begin>a.reference_begin ? parent.reference_begin-a.reference_begin : 0,
                                ob>u ? ob-u : 0);
      const auto end=std::min({a.length,parent.reference_end-a.reference_begin,oe-u});
      for(Length t=start;t<end;) {
        const auto r=a.reference_begin+t;
        const auto q=a.strand==Strand::Forward ? a.query_begin+t : a.query_begin+a.length-1-t;
        auto axis=[&](Position p,Position axis_begin,Position axis_end) {
          struct Spans {std::array<std::pair<Position,Position>,4> values{};std::size_t count{};};
          Spans spans;
          for(unsigned mirror=0;mirror<2;++mirror) {
            const auto position=mirror ? axis_begin+(axis_end-1-p) : p;
            const auto cell=(position-axis_begin)/step;
            for(unsigned previous=0;previous<2;++previous) {
              if(cell<previous)continue;
              const auto low=axis_begin+(cell-previous)*step;
              const auto high=low+std::min(width,axis_end-low);
              const auto span=std::pair{mirror ? axis_begin+(axis_end-high) : low,mirror ? axis_begin+(axis_end-low) : high};
              if(std::find(spans.values.begin(),spans.values.begin()+static_cast<std::ptrdiff_t>(spans.count),span)==
                  spans.values.begin()+static_cast<std::ptrdiff_t>(spans.count))spans.values[spans.count++]=span;
            }
          }
          return spans;
        };
        const auto rs=axis(r,parent.reference_begin,parent.reference_end);
        const auto qs=axis(q,parent.query_begin,parent.query_end);
        for(std::size_t ri=0;ri<rs.count;++ri) for(std::size_t qi=0;qi<qs.count;++qi) {
          const auto [rb,re]=rs.values[ri];const auto [qb,qe]=qs.values[qi];
          DirectionBox b{rb,re,qb,qe};
          boxes.Insert(b,guard);
        }
        const auto rn=step-(r-parent.reference_begin)%step;
        const auto qn=a.strand==Strand::Forward ? step-(q-parent.query_begin)%step :
                                                   (q-parent.query_begin)%step+1;
        const auto rmirror=(parent.reference_end-1-r)%step+1;
        const auto qmirror=a.strand==Strand::Forward ? (parent.query_end-1-q)%step+1 :
                                                       step-(parent.query_end-1-q)%step;
        t+=std::min({rn,qn,rmirror,qmirror,end-t});
      }
    }
    if(guard.options.memory_diagnostic_callback)guard.options.memory_diagnostic_callback(
      "stage=window-scale width="+std::to_string(width)+" windows="+std::to_string(boxes.size())+
      " window_capacity_bytes="+std::to_string(boxes.bytes()));
    Check(guard.options);
    if(width>extent/2) break;
    width*=2;
  }
  return boxes.Finish();
}
std::uint64_t ChainBytes(const DirectionChain& c) {
  std::uint64_t bytes=static_cast<std::uint64_t>(c.nodes.capacity())*sizeof(DirectionChainNode);
  for(const auto& n:c.nodes) bytes+=static_cast<std::uint64_t>(n.source_ids.capacity())*sizeof(std::uint64_t);
  return bytes;
}

// Local best paths generate event hypotheses; every endpoint optimum remains
// in the exported graph. Event coordinates are anchor, not window, bounds.
void Classify(CentromereDirectionResult& out,MemoryGuard& guard,
              std::string_view reference,std::string_view query,std::string_view reverse_query) {
  using Key=std::tuple<DirectionBox,Strand>;
  std::map<Key,InversionEvidence> blocks;
  for(const auto& region:out.regions) {
    Check(guard.options);
    for(std::size_t d=0;d<2;++d) {
      const auto strand=d==0 ? Strand::Forward : Strand::Reverse;
      const auto& chain=region.chains[d];
      if(chain.best==kNoDirectionNode) continue;
      for(const auto end:{chain.best}) {
        // Each local window contributes its preferred path in EACH direction.
        // All other endpoint paths are exported in chains.tsv, not expanded
        // into quadratically many duplicate event hypotheses.
        const auto path=Trace(chain,end);
        if(path.empty()) continue;
        InversionEvidence e;
        e.region_id=region.id; e.strand=strand;
        e.background=d==0 ? Strand::Reverse : Strand::Forward;
        e.support=chain.nodes[end].support;
        const auto first=NodeBox(chain.nodes[path.front()],strand,out.query_length);
        const auto last=NodeBox(chain.nodes[path.back()],strand,out.query_length);
        e.block={first.reference_begin,last.reference_end,
                 std::min(first.query_begin,last.query_begin),std::max(first.query_end,last.query_end)};
        for(const auto i:path)
          e.anchors.insert(e.anchors.end(),chain.nodes[i].source_ids.begin(),chain.nodes[i].source_ids.end());
        std::sort(e.anchors.begin(),e.anchors.end());
        e.anchors.erase(std::unique(e.anchors.begin(),e.anchors.end()),e.anchors.end());
        const auto key=Key{e.block,e.strand};
        auto it=blocks.find(key);
        if(it==blocks.end()) {
          guard.Add(sizeof(InversionEvidence)+64+e.anchors.capacity()*sizeof(std::uint64_t));
          blocks.emplace(key,std::move(e));
        } else if(e.support>it->second.support) it->second=std::move(e);
      }
    }
  }
  const DirectionRegion* root=nullptr;
  for(const auto& r:out.regions) if(r.depth==0 && r.box==out.tasks.front().box) {root=&r;break;}
  for(auto& [key,e]:blocks) {
    (void)key; Check(guard.options);
    const auto opposite=ScoreDirectionChain(out.anchors,e.block,out.reference_length,out.query_length,e.background,
                                            [&]{Check(guard.options);});
    e.opposite_support=opposite.support;
    e.competing=opposite.support>=e.support && opposite.support!=0;
    // Same-direction alternative placements overlapping the same reference
    // block are not allowed to acquire confidence merely by being long.
    for(const auto& [other_key,other]:blocks) {
      (void)other_key;
      if(other.strand!=e.strand || other.support<e.support || other.block==e.block) continue;
      if(other.block.reference_begin>=e.block.reference_end || other.block.reference_end<=e.block.reference_begin) continue;
      if(other.block.query_end<=e.block.query_begin || other.block.query_begin>=e.block.query_end) {
        // Compare equal reference territory. A long alternative sharing only
        // one base must not suppress a short, fully supported local inversion.
        const DirectionBox alternative{e.block.reference_begin,e.block.reference_end,
          other.block.query_begin,other.block.query_end};
        if(ScoreDirectionChain(out.anchors,alternative,out.reference_length,out.query_length,
             e.strand,[&]{Check(guard.options);}).support>=e.support) {e.competing=true;break;}
      }
    }
    const auto bu=e.background==Strand::Forward ? e.block.query_begin : out.query_length-e.block.query_end;
    const auto be=e.background==Strand::Forward ? e.block.query_end : out.query_length-e.block.query_begin;
    const DirectionAnchor *left=nullptr,*right=nullptr;
    auto rank=[&](const DirectionAnchor& a,bool is_left) {
      const auto u=U(a,out.query_length);
      const auto rg=is_left ? e.block.reference_begin-a.reference_begin-a.length : a.reference_begin-e.block.reference_end;
      const auto qg=is_left ? bu-u-a.length : u-be;
      return std::tuple{!(a.full_reference_count==1 && a.full_query_count==1),
                         std::max(rg,qg),rg,qg,std::numeric_limits<Length>::max()-a.length,a.id};
    };
    for(const auto& a:out.anchors) {
      if(a.strand!=e.background) continue;
      const auto u=U(a,out.query_length);
      if(a.reference_begin+a.length<=e.block.reference_begin && u+a.length<=bu &&
         (!left || rank(a,true)<rank(*left,true))) left=&a;
      if(a.reference_begin>=e.block.reference_end && u>=be &&
         (!right || rank(a,false)<rank(*right,false))) right=&a;
    }
    if(left) {
      e.left_flank=left->id; e.left_bounded=true;
      e.reference_left={left->reference_begin+left->length,e.block.reference_begin};
      const auto edge=left->strand==Strand::Forward ? left->query_begin+left->length : left->query_begin;
      const auto boundary=e.background==Strand::Forward ? e.block.query_begin : e.block.query_end;
      e.query_left={std::min(edge,boundary),std::max(edge,boundary)};
    }
    if(right) {
      e.right_flank=right->id; e.right_bounded=true;
      e.reference_right={e.block.reference_end,right->reference_begin};
      const auto edge=right->strand==Strand::Forward ? right->query_begin : right->query_begin+right->length;
      const auto boundary=e.background==Strand::Forward ? e.block.query_end : e.block.query_begin;
      e.query_right={std::min(edge,boundary),std::max(edge,boundary)};
    }
    // Exact anchors can extend across a biological junction by microhomology.
    // An anchor gap therefore need not bound the biological breakpoint tightly.
    // Enumerate equivalent symmetric boundary shifts at the block edges. Both
    // the unchanged and inverted explanations must match on BOTH sequences.
    // This only widens uncertainty; it never adds support or creates anchors.
    if(left && right) {
      const auto oriented=e.background==Strand::Forward ? query : reverse_query;
      auto complement=[](char c) {switch(c){case 'A':return 'T';case 'T':return 'A';case 'C':return 'G';case 'G':return 'C';default:return 'N';}};
      auto equivalent=[&](Position rl,Position rr,Position ql,Position qr) {
        const char a=reference[rl],b=reference[rr];
        return a!='N' && b!='N' && a==complement(b) && oriented[ql]==a && oriented[qr]==b;
      };
      auto widen=[&](Position rl,Position rr,Position ql,Position qr) {
        e.reference_left[0]=std::min(e.reference_left[0],rl);
        e.reference_left[1]=std::max(e.reference_left[1],rl);
        e.reference_right[0]=std::min(e.reference_right[0],rr);
        e.reference_right[1]=std::max(e.reference_right[1],rr);
        const auto l=e.background==Strand::Forward ? ql : out.query_length-ql;
        const auto r=e.background==Strand::Forward ? qr : out.query_length-qr;
        e.query_left[0]=std::min(e.query_left[0],l);e.query_left[1]=std::max(e.query_left[1],l);
        e.query_right[0]=std::min(e.query_right[0],r);e.query_right[1]=std::max(e.query_right[1],r);
      };
      auto rl=e.block.reference_begin,rr=e.block.reference_end,ql=bu,qr=be;
      while(rl>left->reference_begin && rr<right->reference_begin+right->length &&
            ql>U(*left,out.query_length) && qr<U(*right,out.query_length)+right->length &&
            equivalent(rl-1,rr,ql-1,qr)) {--rl;++rr;--ql;++qr;widen(rl,rr,ql,qr);}
      rl=e.block.reference_begin;rr=e.block.reference_end;ql=bu;qr=be;
      while(rr-rl>=2 && qr-ql>=2 && equivalent(rl,rr-1,ql,qr-1)) {
        ++rl;--rr;++ql;--qr;widen(rl,rr,ql,qr);
      }
    }
    const bool unique_left=left && left->full_reference_count==1 && left->full_query_count==1;
    const bool unique_right=right && right->full_reference_count==1 && right->full_query_count==1;
    if(left && right) {
      // The narrow candidate rectangle can hide the unchanged mapping outside
      // its query span. Compare the no-inversion explanation BETWEEN the same
      // flanks before claiming a direction change (essential in repeats).
      const auto lb=U(*left,out.query_length)+left->length;
      const auto rb=U(*right,out.query_length);
      const DirectionBox context{left->reference_begin+left->length,right->reference_begin,
        e.background==Strand::Forward ? lb : out.query_length-rb,
        e.background==Strand::Forward ? rb : out.query_length-lb};
      if(Nonempty(context)) {
        e.context_background_support=ScoreDirectionChain(out.anchors,context,out.reference_length,
          out.query_length,e.background,[&]{Check(guard.options);}).support;
        if(e.context_background_support>=e.support)e.competing=true;
      }
    }
    const auto target=ScoreDirectionChain(out.anchors,e.block,out.reference_length,out.query_length,e.strand);
    const auto actual=target.best==kNoDirectionNode ? 0 : Trace(target,target.best).size();
    e.independent_segments=actual;
    for(const auto i:Trace(target,target.best)) e.longest_match=std::max(e.longest_match,target.nodes[i].length);
    const bool full_reverse=root && e.strand==Strand::Reverse &&
       e.block.reference_begin==0 && e.block.reference_end==out.reference_length &&
       e.block.query_begin==0 && e.block.query_end==out.query_length &&
       root->chains[1].support>root->chains[0].support;
    if(full_reverse && !left && !right) {
      e.classification="whole_reverse_correspondence";e.reason="input_orientation_not_internal_event";
    } else if(e.competing || (left && !unique_left) || (right && !unique_right)) {
      e.classification="repeat_competition";e.reason="opposite_or_flank_context_or_alternative_mapping_or_nonunique_flank";
    } else if(unique_left && unique_right) {
      e.classification="bilateral_inversion_candidate";
      e.reason=actual>=2 ? "two_connections_and_independent_collinear_support" : "two_connections_and_single_exact_block";
    } else if(unique_left || unique_right) {
      e.classification="single_end_inversion_candidate";
      e.reason="one_connection_unresolved";
    } else {
      e.classification="reverse_relationship_unresolved";e.reason="no_resolved_flanking_context";
    }
    // Ordinary background + blocks remain in region/chain evidence, not a
    // misleading list of reverse events. Nested + blocks with - flanks stay.
    if(e.strand==Strand::Forward && !left && !right) continue;
    e.id=out.candidates.size()+1;
    out.candidates.push_back(std::move(e));
  }
  // A block in a terminally supported input orientation is not an inversion just because
  // two inverted repeat matches flank it. A nested event requires an enclosing
  // opposite-direction hypothesis with both junctions supported. Keep the
  // relationship graph, without inferring a unique evolutionary event order.
  // Do NOT choose background by whole-sequence chain length: a large internal
  // inversion can outscore both normal flanks. Require unique mappings at both
  // actual input termini; otherwise leave the background relationship open.
  std::array<bool,2> terminal_left{},terminal_right{};
  for(const auto& a:out.anchors) {
    if(a.full_reference_count!=1 || a.full_query_count!=1)continue;
    const auto d=a.strand==Strand::Forward ? 0U : 1U;
    const auto u=U(a,out.query_length);
    terminal_left[d]=terminal_left[d] || (a.reference_begin==0 && u==0);
    terminal_right[d]=terminal_right[d] || (a.reference_begin+a.length==out.reference_length && u+a.length==out.query_length);
  }
  const bool normal=terminal_left[0] && terminal_right[0],reverse=terminal_left[1] && terminal_right[1];
  if(normal!=reverse) {
    const auto background=normal ? Strand::Forward : Strand::Reverse;
    for(auto& e:out.candidates) {
      if(e.strand!=background || e.classification!="bilateral_inversion_candidate") continue;
      for(const auto& parent:out.candidates) {
        if(parent.strand==e.strand || parent.classification!="bilateral_inversion_candidate" ||
           parent.block==e.block || !Inside(e.block,parent.block)) continue;
        bool contains_flanks=true;
        for(const auto id:{e.left_flank,e.right_flank}) {
          const auto& a=out.anchors.at(static_cast<std::size_t>(id-1));
          contains_flanks=contains_flanks && Inside({a.reference_begin,a.reference_begin+a.length,
            a.query_begin,a.query_begin+a.length},parent.block);
        }
        if(contains_flanks) {guard.Add(sizeof(std::uint64_t)*2);e.enclosing_opposite_candidates.push_back(parent.id);}
      }
      if(e.enclosing_opposite_candidates.empty()) {
        e.classification="reverse_relationship_unresolved";
        e.reason="terminal_input_orientation_without_enclosing_opposite_context";
      }
    }
  }
}
} // namespace
#endif

CentromereDirectionResult AnalyzeCentromereDirections(const SequenceRecord& reference,
    const SequenceRecord& query,const CentromereDirectionOptions& options) {
#if !RAMAG_HAVE_CENTROMERE
  (void)reference;(void)query;(void)options;
  throw DependencyError("bidirectional RaMA backend requires RAMAG_ENABLE_CENTROMERE");
#else
  if(!options.worker_threads || !options.evidence_memory_limit_bytes)
    throw CentromereInputError("direction thread/memory budget must be positive");
  if(options.search_strategy!=DirectionSearchStrategy::RamaFirstHitFocused &&
     options.search_strategy!=DirectionSearchStrategy::LegacyMultiscale)
    throw CentromereInputError("invalid direction search strategy");
  const bool focused=options.search_strategy==DirectionSearchStrategy::RamaFirstHitFocused;
  constexpr Length limit=(INT_MAX-64)/3;
  if(!reference.size() || !query.size() || reference.size()>limit || query.size()>limit-reference.size())
    throw CentromereInputError("direction input empty or exceeds pinned RaMA coordinate range");
  auto complement=[](char c) {switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';case 'N':return 'N';default:throw CentromereInputError("direction input must be normalized A/C/G/T/N");}};
  for(char c:reference.bases) (void)complement(c);
  std::string rc(query.bases.size(),'N');
  for(std::size_t i=0;i<rc.size();++i) rc[rc.size()-1-i]=complement(query.bases[i]);
  const auto start=Clock::now();
  CentromereDirectionResult out;
  out.reference_length=reference.size();out.query_length=query.size();
  MemoryGuard guard{options};
  direction_detail::AnchorLookup known(out.anchors);
  DirectionOriginInterner origins;
  std::uint64_t origin_task=0;
  std::uint64_t callbacks=0,logical_origins=0;
  auto diagnostic=[&](std::string_view stage) {
    if(!options.memory_diagnostic_callback)return;
    std::ostringstream s;
    s<<"stage="<<stage<<" callbacks="<<callbacks<<" anchors="<<out.anchors.size()
      <<" logical_origins="<<logical_origins<<" shared_origin_nodes="<<origins.size()
      <<" anchor_capacity_bytes="<<out.anchors.capacity()*sizeof(DirectionAnchor)
      <<" lookup_capacity_bytes="<<known.bytes()<<" retained_estimate_bytes="<<guard.retained
      <<" peak_estimate_bytes="<<guard.peak<<" tasks="<<out.tasks.size()<<" regions="<<out.regions.size();
    options.memory_diagnostic_callback(s.str());
  };
  std::set<DirectionBox> scheduled;
  auto schedule=[&](const DirectionBox& b,const DirectionTask* parent) {
    if(!Nonempty(b) || (parent && (!Inside(b,parent->box) || b==parent->box))) return;
    if(!scheduled.insert(b).second) return;
    guard.Add(sizeof(DirectionTask)+160);
    out.tasks.push_back({out.tasks.size()+1,parent ? parent->id : 0,parent ? parent->depth+1 : 0,b,{}});
  };
  schedule({0,reference.size(),0,query.size()},nullptr);
  std::size_t begin=0;
  std::vector<direction_detail::CoveredInterval> coverage;
  while(begin<out.tasks.size()) {
    const std::size_t finish=out.tasks.size();
    if(options.progress_callback) options.progress_callback("direction-search",begin,finish);
    std::vector<bool> novel(finish-begin,false);
    // Novelty is defined against the previous breadth-first round, never
    // against updates in this round. Build that immutable union once.
    const auto& previous_coverage=coverage;
    // At most one suffix index resident. Rebuild once per orientation per
    // breadth-first round; no recursive worker pool or nested OpenMP teams.
    for(std::size_t d=0;d<2;++d) {
      Check(options);
      const auto strand=d==0 ? Strand::Forward : Strand::Reverse;
      const std::string_view oriented=d==0 ? std::string_view(query.bases) : std::string_view(rc);
      std::vector<rama_source::SequenceInfo> sequences{
        {reference.bases,reference.name,static_cast<uint_t>(reference.size())},
        {oriented,query.name,static_cast<uint_t>(query.size())}};
      rama_source::AnchorFinder finder(sequences,"",0,false,false,100,[&]{Check(options);});
      std::map<std::string_view,std::pair<uint_t,uint_t>> full_counts;
      for(std::size_t ti=begin;ti<finish;++ti) {
        const auto task=out.tasks[ti];
        // Future task IDs cannot append this task's origin. Drop lookup-only
        // weak entries; live anchor histories keep their own immutable nodes.
        if(origin_task!=task.id) {
          guard.Release(origins.ClearCache()*144);origin_task=task.id;
        }
        const auto& b=task.box;
        finder.visitRegionCandidates({static_cast<uint_t>(b.reference_begin),static_cast<uint_t>(b.reference_end-b.reference_begin),
          static_cast<uint_t>(d==0 ? b.query_begin : query.size()-b.query_end),static_cast<uint_t>(b.query_end-b.query_begin)},
          !focused && options.supplemental_levels,[&](const rama_source::RareCandidate& c) {
            if((++callbacks&1048575ULL)==0)diagnostic("candidate-search");
            const auto& p=c.pair;
            if(!p.match_length || p.second_pos<reference.size()+1) throw AlignmentError("invalid raw RaMA candidate");
            const Length u=p.second_pos-reference.size()-1;
            if(p.first_pos>reference.size() || p.match_length>reference.size()-p.first_pos ||
               u>query.size() || p.match_length>query.size()-u) throw AlignmentError("raw candidate outside input");
            for(Length i=0;i<p.match_length;++i)
              if(reference.bases[p.first_pos+i]=='N' || reference.bases[p.first_pos+i]!=oriented[u+i])
                throw AlignmentError("non-exact raw direction candidate");
            const auto qb=d==0 ? u : query.size()-u-p.match_length;
            if(!Inside({p.first_pos,static_cast<Position>(p.first_pos)+p.match_length,
                        qb,qb+p.match_length},b))
              throw AlignmentError("raw candidate outside its search interval");
            auto index=known.Find(p.first_pos,qb,p.match_length,strand);
            if(index==kNoDirectionNode) {
              const auto pattern=std::string_view(reference.bases).substr(p.first_pos,p.match_length);
              auto count=full_counts.find(pattern);
              if(count==full_counts.end()) {
                count=full_counts.emplace(pattern,finder.countFullOccurrences(p.first_pos,p.match_length)).first;
                guard.Add(80);
              }
              DirectionAnchor a;
              a.id=out.anchors.size()+1;a.reference_begin=p.first_pos;a.query_begin=qb;a.length=p.match_length;
              a.strand=strand;a.source_weight=p.weight;
              a.full_reference_count=count->second.first;a.full_query_count=count->second.second;
              if(out.anchors.size()==out.anchors.capacity()) {
                const auto new_capacity=std::max<std::size_t>(1024,out.anchors.capacity()*2);
                guard.Add(new_capacity*sizeof(DirectionAnchor));
                const auto old=out.anchors.capacity()*sizeof(DirectionAnchor);
                out.anchors.reserve(new_capacity);guard.Release(old);
              }
              index=out.anchors.size();out.anchors.push_back(std::move(a));known.Insert(index,guard);
            }
            auto& a=out.anchors[index];
            const bool was_known=direction_detail::Contains(previous_coverage,a,query.size());
            novel[ti-begin]=!was_known || novel[ti-begin];
            a.source_weight=std::max(a.source_weight,p.weight);
            const DirectionOrigin origin{task.id,task.depth,c.level,c.reference_count,c.query_count,c.first_hit_level};
            if(a.origins.empty() || a.origins.back().task_id!=origin.task_id || a.origins.back().rare_level!=origin.rare_level ||
               a.origins.back().reference_interval_count!=origin.reference_interval_count ||
               a.origins.back().query_interval_count!=origin.query_interval_count) {
              if(origins.Append(a.origins,origin))guard.Add(sizeof(DirectionOriginNode)+144);
              ++logical_origins;
            }
          });
      }
      diagnostic(d==0 ? "forward-search-complete" : "reverse-search-complete");
      guard.Release(full_counts.size()*80);
    }
    guard.Release(origins.ClearCache()*144);
    auto updated_coverage=direction_detail::Coverage(out.anchors,query.size(),guard);
    guard.Release(coverage.capacity()*sizeof(direction_detail::CoveredInterval));
    coverage.swap(updated_coverage);
    std::vector<direction_detail::CoveredInterval>().swap(updated_coverage);
    diagnostic("mapping-union-complete intervals="+std::to_string(coverage.size()));
    for(std::size_t ti=begin;ti<finish;++ti) {
      Check(options);
      const auto task=out.tasks[ti];
      std::vector<direction_detail::CompactBox> boxes;
      if(focused) {
        // One corresponding region per search task. There is no coordinate
        // grid, and no minimum window size or hidden candidate truncation.
        guard.Add(sizeof(direction_detail::CompactBox));
        boxes.push_back(direction_detail::CompactBox(task.box));
      } else boxes=Windows(coverage,task.box,query.size(),guard);
      diagnostic("windows-complete count="+std::to_string(boxes.size()));
      const auto before=out.tasks.size();
      const bool expand=options.recurse && (task.depth==0 || novel[ti-begin]);
      for(const auto& compact_box:boxes) {
        const auto box=compact_box.Box();
        DirectionRegion region;
        region.id=out.regions.size()+1;region.task_id=task.id;region.depth=task.depth;region.box=box;
        region.scale=std::max(box.reference_end-box.reference_begin,box.query_end-box.query_begin);
        for(std::size_t d=0;d<2;++d) {
          const auto strand=d==0 ? Strand::Forward : Strand::Reverse;
          region.chains[d]=ScoreDirectionChain(out.anchors,box,reference.size(),query.size(),strand,[&]{Check(options);});
          guard.Add(ChainBytes(region.chains[d]));
          if(!expand || focused) continue;
          const auto& chain=region.chains[d];
          // Recursive proposals come from each window's best path in BOTH
          // directions. All endpoint alternatives are still retained in the
          // graph; recursively expanding every traceback edge creates the
          // same tiny rectangle combinatorially, without new mappings.
          for(const auto ni:Trace(chain,chain.best)) {
            const auto& n=chain.nodes[ni];
            schedule(NodeBox(n,strand,query.size()),&task); // inspect the opposite orientation even inside a covered match
            const auto ub=d==0 ? box.query_begin : query.size()-box.query_end;
            const auto ue=d==0 ? box.query_end : query.size()-box.query_begin;
            auto gap=[&](Position rb,Position re,Position qb,Position qe) {
              if(rb>=re || qb>=qe)return;
              schedule({rb,re,d==0 ? qb : query.size()-qe,d==0 ? qe : query.size()-qb},&task);
            };
            if(n.predecessor!=kNoDirectionNode) {
              const auto& p=chain.nodes[n.predecessor];
              gap(p.reference_begin+p.length,n.reference_begin,p.oriented_query_begin+p.length,n.oriented_query_begin);
            } else gap(box.reference_begin,n.reference_begin,ub,n.oriented_query_begin);
            gap(n.reference_begin+n.length,box.reference_end,n.oriented_query_begin+n.length,ue);
          }
        }
        const auto f=region.chains[0].support,r=region.chains[1].support;
        const double sum=static_cast<double>(f)+static_cast<double>(r);
        region.advantage=sum ? (static_cast<double>(f)-static_cast<double>(r))/sum : 0;
        if(focused && expand) {
          // First-hit RaMA gap recursion in BOTH orientations. A path emits
          // each adjacent gap once, plus its two terminal gaps, never a tail
          // rectangle for every node. The opposite direction remains evidence
          // and is searched again within every proposed child.
          for(std::size_t d=0;d<2;++d) {
            const auto strand=d==0 ? Strand::Forward : Strand::Reverse;
            const auto& chain=region.chains[d];
            // A global length-optimal path may cross short chance matches.
            // Keep its local coherent branches as well. A join whose original
            // RaMA discordance cost exhausts either adjacent run's source
            // weight starts a new block; equal diagonals are always coherent.
            // These blocks propose regions, not a new direction score.
            std::vector<std::size_t> first(chain.nodes.size());
            std::vector<bool> extended(chain.nodes.size(),false);
            for(std::size_t i=0;i<chain.nodes.size();++i) {
              first[i]=i;const auto& n=chain.nodes[i];
              if(n.predecessor==kNoDirectionNode)continue;
              const auto& p=chain.nodes[n.predecessor];
              const auto rg=n.reference_begin-p.reference_begin,qg=n.oriented_query_begin-p.oriented_query_begin;
              const auto delta=rg>=qg ? rg-qg : qg-rg;
              const auto cost=2.0*std::log2(static_cast<double>(delta)+1.0);
              if(!delta || cost<std::min(n.weight,p.weight)) {
                first[i]=first[n.predecessor];extended[n.predecessor]=true;
              }
            }
            for(std::size_t i=0;i<chain.nodes.size();++i) {
              if(extended[i])continue;
              const auto a=NodeBox(chain.nodes[first[i]],strand,query.size());
              const auto z=NodeBox(chain.nodes[i],strand,query.size());
              DirectionBox block{a.reference_begin,z.reference_end,
                std::min(a.query_begin,z.query_begin),std::max(a.query_end,z.query_end)};
              schedule(block,&task);
              // Junction microhomology can make an exact run extend into
              // its opposite-direction flank. Also inspect the rectangle
              // between those actual paired flanks; keep the untrimmed block.
              const auto background=d==0 ? Strand::Reverse : Strand::Forward;
              const auto bu=background==Strand::Forward ? block.query_begin : query.size()-block.query_end;
              const auto be=background==Strand::Forward ? block.query_end : query.size()-block.query_begin;
              const DirectionChainNode *left=nullptr,*right=nullptr;
              for(const auto& p:region.chains[1-d].nodes) {
                const auto re=p.reference_begin+p.length,qe=p.oriented_query_begin+p.length;
                auto better=[](const DirectionChainNode& candidate,const DirectionChainNode* old) {
                  return !old || candidate.length>old->length || (candidate.length==old->length &&
                    std::tie(candidate.reference_begin,candidate.oriented_query_begin)<std::tie(old->reference_begin,old->oriented_query_begin));
                };
                if(p.reference_begin<=block.reference_begin && block.reference_begin<re && re<block.reference_end &&
                   p.oriented_query_begin<=bu && bu<qe && qe<be && better(p,left))left=&p;
                if(block.reference_begin<p.reference_begin && p.reference_begin<block.reference_end && block.reference_end<=re &&
                   bu<p.oriented_query_begin && p.oriented_query_begin<be && be<=qe && better(p,right))right=&p;
              }
              auto trimmed=block;auto low=bu,high=be;
              if(left) {trimmed.reference_begin=left->reference_begin+left->length;low=left->oriented_query_begin+left->length;}
              if(right) {trimmed.reference_end=right->reference_begin;high=right->oriented_query_begin;}
              if(low<high) {
                trimmed.query_begin=background==Strand::Forward ? low : query.size()-high;
                trimmed.query_end=background==Strand::Forward ? high : query.size()-low;
                schedule(trimmed,&task);
              }
            }
            const auto path=Trace(chain,chain.best);
            if(path.empty())continue;
            auto rb=box.reference_begin;
            auto ub=d==0 ? box.query_begin : query.size()-box.query_end;
            const auto ue=d==0 ? box.query_end : query.size()-box.query_begin;
            auto gap=[&](Position re,Position qe) {
              if(rb<re && ub<qe)
                schedule({rb,re,d==0 ? ub : query.size()-qe,d==0 ? qe : query.size()-ub},&task);
            };
            for(const auto ni:path) {
              const auto& n=chain.nodes[ni];
              gap(n.reference_begin,n.oriented_query_begin);
              rb=n.reference_begin+n.length;ub=n.oriented_query_begin+n.length;
              // A covered exact run is only revisited when the other
              // orientation has actual paired-column support inside it.
              const auto b=NodeBox(n,strand,query.size());
              const auto opposite=d==0 ? Strand::Reverse : Strand::Forward;
              const auto ob=opposite==Strand::Forward ? b.query_begin : query.size()-b.query_end;
              const auto oe=opposite==Strand::Forward ? b.query_end : query.size()-b.query_begin;
              for(const auto& p:region.chains[1-d].nodes) {
                if(p.reference_begin>=b.reference_end)break;
                if(p.oriented_query_begin>=oe)continue;
                const auto start_offset=std::max(b.reference_begin>p.reference_begin ? b.reference_begin-p.reference_begin : 0,
                  ob>p.oriented_query_begin ? ob-p.oriented_query_begin : 0);
                const auto stop_offset=std::min({p.length,b.reference_end-p.reference_begin,oe-p.oriented_query_begin});
                if(start_offset<stop_offset) {schedule(b,&task);break;}
              }
            }
            gap(box.reference_end,ue);
          }
        } else if(expand && (f || r) && region.scale<=80) schedule(box,&task);
        guard.Add(2*sizeof(DirectionRegion));out.regions.push_back(std::move(region));
      }
      out.tasks[ti].stop_reason=!options.recurse ? "recursion_disabled" :
        !expand ? "no_novel_mapping" : out.tasks.size()==before ? "no_new_smaller_region" : "children_scheduled";
      if(options.progress_callback) options.progress_callback("centromere-directions",ti+1,out.tasks.size());
      guard.Release(boxes.capacity()*sizeof(direction_detail::CompactBox));
    }
    begin=finish;
  }
  guard.Release(coverage.capacity()*sizeof(direction_detail::CoveredInterval));
  std::vector<direction_detail::CoveredInterval>().swap(coverage);
  if(options.progress_callback) options.progress_callback("direction-classify",0,out.regions.size());
  Classify(out,guard,reference.bases,query.bases,rc);
  diagnostic("complete");
  auto& s=out.statistics;
  s.reference_contigs=s.query_contigs=1;s.reference_bases=reference.size();s.query_bases=query.size();
  s.actual_seed_route=std::string(DirectionStrategyName(options.search_strategy));
  s.chaining_route="effective-collinear-length-v1";
  s.seed_requested_threads=options.worker_threads;s.seed_scheduled_threads=s.seed_worker_threads=1;
  s.seed_parallel_route="serial-direction-rounds-no-nested-workers";
  s.total_seconds=std::chrono::duration<double>(Clock::now()-start).count();s.seed_seconds=s.total_seconds;
  return out;
#endif
}
} // namespace ramag

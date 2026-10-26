#pragma once
#include "centromere_orientation_connections.hpp"
#include <memory_resource>
#include <vector>

namespace ramag::detail {
// One frozen graph per original long gap. Its vertices are actual exact-query
// unions and the two surviving old endpoints, never inferred interval edges.
// This is direction evidence, not a chain of compatible reference-copy mappings.
template<class EmitPiece,class EmitGap,class Check>
void VisitSparseOrientationPivots(std::span<const OrientationEndpoint> endpoints,
    std::span<const OrientationCoverage> labels,std::string_view query,Length limit,
    std::span<const CentromereQueryInterval> known_forward,
    std::span<const CentromereQueryInterval> known_reverse,
    std::span<const CentromereQueryInterval> forward,
    std::span<const CentromereQueryInterval> reverse,
    std::span<const CentromereQueryInterval> occupied,std::pmr::memory_resource* memory,
    EmitPiece emit_piece,EmitGap emit_gap,Check check) {
  if(!limit || !memory)throw AlignmentError("sparse-pivot recovery requires a positive limit and workspace");
  auto validate=[&](auto rows) {
    for(std::size_t i=0;i<rows.size();++i) {
      check();const auto& row=rows[i];
      if(row.begin>=row.end || row.end>query.size() || (i && rows[i-1].end>row.begin))
        throw AlignmentError("invalid sparse-pivot interval index");
    }
  };
  validate(endpoints);validate(labels);validate(known_forward);validate(known_reverse);
  validate(forward);validate(reverse);validate(occupied);
  for(const auto& row:labels)if(row.mask>7U)throw AlignmentError("invalid sparse-pivot window mask");
  auto contains=[](auto rows,Position p) {
    const auto it=std::lower_bound(rows.begin(),rows.end(),p,[](const auto& row,auto x){return row.end<=x;});
    return it!=rows.end() && it->begin<=p;
  };
  using Interval=CentromereQueryInterval;
  std::pmr::vector<Interval> vertices(memory);
  struct Edge {bool pass{};std::string_view reason;};
  std::pmr::vector<Edge> edges(memory);
  for(std::size_t i=1;i<endpoints.size();++i) {
    check();const auto& left=endpoints[i-1];const auto& right=endpoints[i];
    if(left.strand!=right.strand || right.begin-left.end<=limit)continue;
    if(left.strand!=Strand::Forward && left.strand!=Strand::Reverse)
      throw AlignmentError("invalid sparse-pivot endpoint strand");
    const Position begin=left.end,end=right.begin;
    const auto own=left.strand==Strand::Forward?forward:reverse;
    const auto opposite=left.strand==Strand::Forward?reverse:forward;
    const auto known=left.strand==Strand::Forward?known_forward:known_reverse;
    const unsigned bit=left.strand==Strand::Forward?1U:2U;
    auto decision=[&](Position b,Position e,bool accepted,std::string_view reason) {
      emit_gap(CentromereSparseOrientationGap{b,e,begin,end,left.strand,accepted,std::string(reason)});
    };
    if(!contains(known,begin-1) || !contains(own,begin-1) ||
       !contains(known,end) || !contains(own,end)) {
      decision(begin,end,false,"original_pivot_not_supported");continue;
    }
    auto barrier=[&](Position b,Position e)->std::string_view {
      if(b==e)return {};
      auto label=std::lower_bound(labels.begin(),labels.end(),b,[](const auto& row,auto p){return row.end<=p;});
      for(;label!=labels.end() && label->begin<e;++label) {
        check();if(label->mask & ~bit)return "opposing_or_competitive_window";
      }
      const auto other=std::lower_bound(opposite.begin(),opposite.end(),b,[](const auto& row,auto p){return row.end<=p;});
      if(other!=opposite.end() && other->begin<e)return "opposite_exact_evidence";
      if(query.substr(b,e-b).find('N')!=std::string_view::npos)return "contains_N";
      return {};
    };
    vertices.clear();edges.clear();vertices.push_back({begin,begin});
    auto root=std::lower_bound(own.begin(),own.end(),begin,[](const auto& row,auto p){return row.end<=p;});
    for(;root!=own.end() && root->begin<end;++root) {
      check();const auto b=std::max(begin,root->begin),e=std::min(end,root->end);
      const auto reason=barrier(b,e);
      // Never manufacture shorter, newly eligible pivots by removing a vetoed
      // base from a root segment. The normalized query segment is all-or-nothing.
      if(!reason.empty()){decision(b,e,false,reason);continue;}
      vertices.push_back({b,e});
    }
    if(vertices.size()==1){decision(begin,end,false,"no_valid_internal_root");continue;}
    vertices.push_back({end,end});
    for(std::size_t j=1;j<vertices.size();++j) {
      check();const auto b=vertices[j-1].end,e=vertices[j].begin;
      const auto reason=e-b>limit?std::string_view{"gap_too_long"}:barrier(b,e);
      edges.push_back({reason.empty(),reason});
    }
    // A linear graph needs only two reachability frontiers, not path copies or
    // recursive rescans. An isolated middle component cannot claim direction.
    std::size_t l=0,r=vertices.size()-1;
    while(l<edges.size() && edges[l].pass){check();++l;}
    while(r && edges[r-1].pass){check();--r;}
    auto emit_unoccupied=[&](Position b,Position e,bool inferred) {
      auto used=std::lower_bound(occupied.begin(),occupied.end(),b,[](const auto& row,auto p){return row.end<=p;});
      Position cursor=b;
      for(;used!=occupied.end() && used->begin<e;++used) {
        check();
        if(cursor<used->begin)emit_piece(CentromereSparseOrientationInterval{
          0,cursor,std::min(e,used->begin),begin,end,left.strand,inferred});
        cursor=std::max(cursor,used->end);
        if(cursor>=e)return;
      }
      if(cursor<e)emit_piece(CentromereSparseOrientationInterval{0,cursor,e,begin,end,left.strand,inferred});
    };
    for(std::size_t j=1;j<vertices.size();++j) {
      check();const auto b=vertices[j-1].end,e=vertices[j].begin;
      const bool connected=edges[j-1].pass && (j<=l || j-1>=r);
      if(b<e)decision(b,e,connected,connected?"sparse_root_bridge":
        edges[j-1].pass?"isolated_internal_component":edges[j-1].reason);
      // Test the original full edge before subtracting old predictions. Occupied
      // inferred ends cannot turn an over-limit edge into an admissible one.
      if(connected)emit_unoccupied(b,e,true);
      if(j+1<vertices.size() && (j<=l || j>=r))emit_unoccupied(vertices[j].begin,vertices[j].end,false);
    }
  }
}
} // namespace ramag::detail

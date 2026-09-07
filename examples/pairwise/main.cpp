#include <ramag/pairwise_core.hpp>
#include <iostream>

int main() {
    std::vector<ramag::SequenceRecord> reference{{0,"ref","ref",std::string(80,'A')}};
    std::vector<ramag::SequenceRecord> query{{0,"query","query",std::string(80,'A')}};
    std::vector<ramag::Seed> anchors{{0,0,0,0,80,ramag::Strand::Forward}};
    const auto result=ramag::AlignPairwiseCore(reference,query,anchors);
    if(result.records.size()!=1 || result.records.front().record.score!=720) return 1;
    std::cout << "records=" << result.records.size() << '\n';
}

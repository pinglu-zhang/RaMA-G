#include <ramag/pairwise_core.hpp>
#include <ramag/centromere.hpp>
#include <ramag/centromere_direction.hpp>
#include <ramag/centromere_blocks.hpp>
#include <ramag/centromere_directed_anchors.hpp>
#include <ramag/centromere_copy_anchors.hpp>
#include <ramag/centromere_stitch.hpp>
#include <ramag/centromere_stitch_formats.hpp>
#include <iostream>

int main() {
  std::vector<ramag::SequenceRecord> r{{0, "r", "r", std::string(80, 'A')}};
  std::vector<ramag::SequenceRecord> q{{0, "q", "q", std::string(80, 'A')}};
  std::vector<ramag::Seed> seeds{{0, 0, 0, 0, 80, ramag::Strand::Forward}};
  const auto aligned = ramag::AlignPairwiseCore(r, q, seeds);
  if (aligned.records.size() != 1 || aligned.records.front().record.score != 720) return 1;
  if (ramag::CentromereBackendAvailable() != static_cast<bool>(EXPECT_BACKEND)) return 2;
  ramag::CentromereStitchOptions stitch;
  if (stitch.max_dp_cells != 4000000) return 3;
  std::cout << "installed headers, core score and backend contract passed\n";
}

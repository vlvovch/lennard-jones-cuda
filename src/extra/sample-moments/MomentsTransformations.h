/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#ifndef SAMPLEMOMENTS_MOMENTSTRANSFORMATIONS_H
#define SAMPLEMOMENTS_MOMENTSTRANSFORMATIONS_H

#include <vector>
#include <map>
#include <algorithm>
#include <cassert>

namespace SampleMoments {

  /// Returns all partitions of a set {1,2,...,n} via a recursive procedure
  /// Each element contains a vector of blocks, each block is a vector of indices
  /// See https://en.wikipedia.org/wiki/Partition_of_a_set
  /// \param n  The dimension of the set
  static std::vector<std::vector<std::vector<int>>> PartitionsOfSet(int n) {
    if (n == 1)
      return {{{0}}};

    auto subpartitions = PartitionsOfSet(n - 1);
    auto ret = std::vector<std::vector<std::vector<int>>>();

    for (auto &el : subpartitions) {
      auto temp_element = el;
      temp_element.push_back({n - 1});
      ret.push_back(temp_element);
      temp_element = el;
      for (int i = 0; i < el.size(); ++i) {
        temp_element[i].push_back(n - 1);
        ret.push_back(temp_element);
        temp_element[i] = el[i];
      }
    }

    return ret;
  }

  /// Express an arbitrary joint cumulant \kappa(X_i0,X_i1,X_i2,...) in terms of joint (central) moments.
  /// Uses the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
  /// \param indices            A vector of 0-based indices of the joint cumulant
  /// \param dimensions         The total number of distinct random variable. 
  ///                           If set to -1, its minimum possible value is inferred from the contents of the provided cumulant \p indices
  /// \param ordinary_moments   Whether to express the cumulants in terms ordinary moments (true) or central moments (false). The default is central moments
  ///
  /// \return                   Returns the joint cumulant as a sum of products of various (central) moments.
  ///                           The return value is a map. 
  ///                           Each key-value pair correspond to a single term.
  ///                           The key corresponds to a vector of (central) moments. Each element of this vector is a vector of the moments' indices.
  ///                           The value is a numeral factor in front of each term. 
  static std::map<std::vector<std::vector<int>>, int64_t> JointCumulantToCentralMoments(
          const std::vector<int> &indices,
          int dimensions = -1,
          bool ordinary_moments = false
  ) {
    assert(indices.size() > 1);

    if (dimensions == -1) {
      for (int ind : indices) {
        assert(ind >= 0);
        dimensions = std::max(dimensions, ind);
      }
      dimensions++;
    }

    // Use the multivariate version of the Faa di Bruno's formula https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula
    auto ret = std::map<std::vector<std::vector<int>>, int64_t>();

    std::vector<int64_t> factorials(indices.size(), 1);
    for (int k = 1; k < factorials.size(); ++k)
      factorials[k] = k * factorials[k - 1];

    auto partitions_FdB = PartitionsOfSet(indices.size());
    for (const auto &partition : partitions_FdB) {
      int64_t multiplier = 1;
      int num_blocks = partition.size();
      if (num_blocks % 2 == 0)
        multiplier = -1;
      multiplier *= factorials[num_blocks - 1];

      std::vector<std::vector<int>> mults;
      bool zero_term = false;
      for (const auto &block : partition) {
        if (block.size() == 1) {
          zero_term = true;
        }

        auto inds = std::vector<int>(dimensions, 0);
        for (int ind : block)
          inds[indices[ind]]++;

        mults.push_back(inds);
      }
      if (zero_term && !ordinary_moments)
        continue;
      sort(mults.begin(), mults.end());
      ret[mults] += multiplier;
    }

    return ret;
  }

  /// Express an arbitrary joint cumulant in terms of the joint moments
  /// Calls JointCumulantToCentralMoments()
  static std::map<std::vector<std::vector<int>>, int64_t> JointCumulantToMoments(
          const std::vector<int> &indices,
          int dimensions = -1
  ) {
    return JointCumulantToCentralMoments(indices, dimensions, true);
  }

} // namespace SampleMoments

#endif
/*
 * sample-moments library
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>
 */
#ifndef SAMPLEMOMENTS_NUMBERSTATISTICS_H
#define SAMPLEMOMENTS_NUMBERSTATISTICS_H

#include "MomentsTransformations.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <cstdint>

namespace SampleMoments {
  /**
   * \brief A class implementing calculation of various sample statistics, and their standard errors, for a single variable.
   *        For example this includes (higher-order) moments, central moments, cumulants, and ratios of such quantities.
   *
   *
   */
  class NumberStatistics {
    /// The maximum order of moments stored
    int32_t m_Nmax;

    /// The cumulative sums of N^k values collected from all observations
    std::vector<double> m_MomentSums;

    /// The sample means of the moments, <N^k>. Populated by CalculateMoments()
    std::vector<double> m_Moments;

    /// The sample means of the central moments, <(N-<N>)^k>. Populated by CalculateMoments()
    std::vector<double> m_CentralMoments;

    /// Shifts all moments by a consant value
    /// Can be useful to avoid large round-off errors if the expected mean is known
    /// Does not affect the central moments
    double m_MeanShift;

    /// The total number of observations
    int64_t m_NumberOfObservations;

    /// Same as m_NumberOfObservations
    double m_nE;

    /// Whether the moments have been computed from the sums
    bool m_MomentsComputed;

    /// The pre-computed binomial coefficients
    std::vector<std::vector<int64_t> > m_BinomialCoefficients;
  public:
    /**
     * \brief Construct a new NumberStatistics object.
     * 
     * \param nmax           The maximum order of moments to be stored
     * \param observations   The vector of
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    NumberStatistics(int nmax = 16) {
      Reset(nmax);
    }

    /**
     * \brief Construct a new NumberStatistics object.
     *
     * \param observations  A vector of observations
     * \param nmax          The maximum order of moments to be stored
     *
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    template<typename T = double>
    NumberStatistics(const std::vector<T> &observations = std::vector<T>(), int nmax = 16) {
      Reset(nmax);
      AddObservations(observations);
    }

    /// Destructor. Does nothing
    virtual ~NumberStatistics() {}

    /**
     * \brief Resets the observations to zero.
     * 
     * \param nmax   The maximum order of moments to be stored
     * 
     * For the error estimation nmax must be at least twice larger than the order of the moment to be estimated
     */
    virtual void Reset(int nmax) {
      assert(nmax >= 1);
      m_Nmax = nmax;
      m_MomentSums = std::vector<double>(nmax + 1, 0.);
      m_Moments = std::vector<double>(nmax + 1, 0.);
      m_CentralMoments = std::vector<double>(nmax + 1, 0.);
      m_NumberOfObservations = 0;
      m_nE = 0.;
      m_MeanShift = 0.;

      // Pre-compute the binomial coefficients
      CalculateBinomialCoefficients();

      m_MomentsComputed = false;
    }

    /**
     * \brief Clears everything
     * 
     */
    virtual void Clear() {
      Reset(m_Nmax);
    }

    /// Total number of observations
    int64_t GetNumberOfObservations() const { return m_NumberOfObservations; }

    /// Whether the shift of the mean value is applied
    bool IsMeanShifted() {
      return (m_MeanShift != 0.0);
    }

    /// Adds an observation with an integer value of N. This value is cast to double-precision floating-point number
    template<typename T>
    void AddObservation(const T& N) {
      AddObservation(static_cast<double>(N));
    }

    /// Adds an observation with an real value of N
    void AddObservation(double N) {
      UpdateMomentSumsWithObservation(N, m_MomentSums);

      m_NumberOfObservations++;

      m_MomentsComputed = false;
    }

    /// Adds several integer-valued observations
    /// \param observations A vector of observations
    template<typename T>
    void AddObservations(const std::vector<T> &observations) {
      std::vector<double> observations_to_double;
      observations_to_double.reserve(observations.size());
      for (const T &N : observations)
        observations_to_double.push_back(static_cast<double>(N));
      AddObservations(observations_to_double);
    }

    /// Adds several real-valued observations
    /// \param observations A vector of observations
    void AddObservations(const std::vector<double> &observations) {
      // Compute moments_sum from the whole observations vector to minimize round-off errors
      auto moment_sums = std::vector<double>(m_Nmax, 0.);

      for (const double &N : observations) {
        UpdateMomentSumsWithObservation(N, moment_sums);
      }

      // Now add to the whole sums
      for (int i = 0; i < m_Nmax + 1; ++i) {
        m_MomentSums[i] += moment_sums[i];
      }

      m_NumberOfObservations += static_cast<int64_t>(observations.size());

      m_MomentsComputed = false;
    }


    /// Returns the r-th moment <N^k>
    double GetMoment(int r) {
      assert(r <= m_Nmax);
      CalculateMoments();
      return m_Moments[r];
    }


    /// Returns the sample covariance cov(<N^r>,<N^q>)
    /// In accordance with M. Kendall and A. Stuart, ``The advanced theory of statistics''
    double GetMomentsSampleCovariance(int r, int q) {
      assert(r + q <= m_Nmax);
      CalculateMoments();
      return (m_Moments[r + q] - m_Moments[r] * m_Moments[q]) / m_nE;
    }

    /// Returns the error estimate for the r-th moment <N^k>
    double GetMomentError(int r) {
      return std::sqrt(GetMomentsSampleCovariance(r, r));
    }

    /// Returns the ratio of moments <N^r> / <N^q>
    double GetMomentRatio(int r, int q) { return GetMoment(r) / GetMoment(q); }

    /// Returns the error estimate for the ratio of moments <N^r> / <N^q>
    double GetMomentRatioError(int r, int q) {
      double c1 = GetMoment(r);
      double c2 = GetMoment(q);
      double c1Dev = GetMomentsSampleCovariance(r, r);
      double c2Dev = GetMomentsSampleCovariance(q, q);
      double c1c2cov = GetMomentsSampleCovariance(r, q);

      return std::abs(c1 / c2) * std::sqrt(c1Dev / c1 / c1 + c2Dev / c2 / c2 - 2. * c1c2cov / c1 / c2);
    }

    /// Returns the r-th central moment, \mu_r = <(N-<N>)^r>
    double GetCentralMoment(int r) {
      assert(r <= m_Nmax);
      CalculateMoments();
      return m_CentralMoments[r];
    }

    /// Returns the sample covariance cov(<N^r>,<N^q>)
    /// In accordance with M. Kendall and A. Stuart, ``The advanced theory of statistics''
    double GetCentralMomentsSampleCovariance(int r, int q) {
      CalculateMoments();

      assert(m_Nmax >= 2);
      assert(r + q <= m_Nmax);
      if (r <= 1 || q <= 1)
        return 0.0;

      double ret = 0.0;

      ret += m_CentralMoments[r + q];
      ret -= m_CentralMoments[r] * m_CentralMoments[q];
      ret += r * q * m_CentralMoments[2] * m_CentralMoments[r - 1] * m_CentralMoments[q - 1];
      ret -= r * m_CentralMoments[r - 1] * m_CentralMoments[q + 1];
      ret -= q * m_CentralMoments[r + 1] * m_CentralMoments[q - 1];

      return ret / m_nE;
    }

    /// Returns the error estimate for the r-th central moment
    double GetCentralMomentError(int r) {
      return std::sqrt(GetCentralMomentsSampleCovariance(r, r));
    }

    /// Returns the ratio of central moments \mu_r / \mu_q
    double GetCentralMomentRatio(int r, int q) { return GetCentralMoment(r) / GetCentralMoment(q); }

    /// Returns the error estimate for the ratio of central moments \mu_r / \mu_q
    double GetCentralMomentRatioError(int r, int q) {
      double c1 = GetCentralMoment(r);
      double c2 = GetCentralMoment(q);
      double c1Dev = GetCentralMomentsSampleCovariance(r, r);
      double c2Dev = GetCentralMomentsSampleCovariance(q, q);
      double c1c2cov = GetCentralMomentsSampleCovariance(r, q);

      return std::abs(c1 / c2) * std::sqrt(c1Dev / c1 / c1 + c2Dev / c2 / c2 - 2. * c1c2cov / c1 / c2);
    }

    /// Returns the r-th order cumulant \kappa_r
    double GetCumulant(int r) {
      return CalculateCumulantFromMoments(r);
    }

    /// Returns the error estimate for the r-th order cumulant \kappa_r
    double GetCumulantError(int r) {
      return CalculateCumulantErrorFromMoments(r);
    }

    /// Returns the cumulant ratio \kappa_r / \kappa_q
    double GetCumulantRatio(int r, int q) { return GetCumulant(r) / GetCumulant(q); }

    /// Returns the error estimate for the cumulant ratio \kappa_r / \kappa_q
    double GetCumulantRatioError(int r, int q) { return CalculateCumulantRatioErrorFromMoments(r, q); }


    /// Returns the r-th cumulant scaled by the mean \kappa_r / <N>
    double GetScaledCumulant(int r) { return GetCumulant(r) / GetMean(); }

    /// Returns the error estimate for r-th cumulant scaled by the mean \kappa_r / <N>
    double GetScaledCumulantError(int r) {
      if (r == 1)
        return 0.0;

      double c1 = GetCentralMoment(r);
      double c2 = GetMean();
      double c1Dev = GetCentralMomentsSampleCovariance(r, r);
      double c2Dev = GetMeanError() * GetMeanError();

      CalculateMoments();

      double c1c2cov = (m_CentralMoments[r + 1] - r * m_CentralMoments[2] * m_CentralMoments[r - 1]) / m_nE;

      return std::abs(c1 / c2) * std::sqrt(c1Dev / c1 / c1 + c2Dev / c2 / c2 - 2. * c1c2cov / c1 / c2);
    }

    /// Returns the sample mean
    double GetMean() { return GetMoment(1); }

    /// Returns the error estimate of the mean
    double GetMeanError() { return GetMomentError(1); }

    /// Returns the variance of the distribution
    double GetVariance() { return GetCentralMoment(2); }

    /// Returns the error estimate for the variance of the distribution
    double GetVarianceError() { return GetCentralMomentError(2); }

    double GetScaledVariance() { return GetScaledCumulant(2); }

    double GetScaledVarianceError() { return GetScaledCumulantError(2); }


    double GetSkewness() { return GetCumulantRatio(3, 2); }

    double GetSkewnessError() { return GetCumulantRatioError(3, 2); }

    double GetKurtosis() { return GetCumulantRatio(4, 2); }

    double GetKurtosisError() { return GetCumulantRatioError(4, 2); }

    double GetC5C2() { return GetCumulantRatio(5, 2); }

    double GetC5C2Error() { return GetCumulantRatioError(5, 2); }

    double GetC6C2() { return GetCumulantRatio(6, 2); }

    double GetC6C2Error() { return GetCumulantRatioError(6, 2); }

    /// Set the shift of means to use when accumulating statistics
    /// Can be useful to avoid round-off error when dealing with higher-order central moments and cumulants
    void SetMeanShift(double meanshift) {
      if (m_NumberOfObservations == 0) {
        m_MeanShift = meanshift;
        return;
      }

      auto shifted_moments = m_Moments;
      for (int i = 0; i < m_Nmax + 1; ++i)
        shifted_moments[i] = m_MomentSums[i] / m_NumberOfObservations;

      double delta = m_MeanShift - meanshift;
      m_MeanShift = meanshift;

      auto new_shifted_moments = shifted_moments;
      new_shifted_moments[0] = 1.;
      for (int n = 1; n <= m_Nmax; ++n) {
        new_shifted_moments[n] = 0.;
        double Nshift = 1.;
        for (int k = 0; k <= n; ++k) {
          new_shifted_moments[n] += BinomialCoefficient(n, k) * shifted_moments[n - k] * Nshift;
          Nshift *= delta;
        }
      }

      for (int i = 0; i < m_Nmax + 1; ++i) {
        m_MomentSums[i] = new_shifted_moments[i] * m_NumberOfObservations;
      }

      m_MomentsComputed = false;
      CalculateMoments();
    }

  protected:
    /// Pre-compute the needed binomial coefficients
    void CalculateBinomialCoefficients() {
      m_BinomialCoefficients = std::vector<std::vector<int64_t> >(m_Nmax + 1);

      m_BinomialCoefficients[0] = std::vector<int64_t>(1, 1);
      m_BinomialCoefficients[1] = std::vector<int64_t>(2, 1);

      for (int n = 2; n < m_Nmax + 1; ++n) {
        m_BinomialCoefficients[n] = std::vector<int64_t>(n + 1, 1);
        m_BinomialCoefficients[n][0] = m_BinomialCoefficients[n][n] = 1;
        for (int mm = 1; mm < n; ++mm)
          m_BinomialCoefficients[n][mm] = m_BinomialCoefficients[n - 1][mm] + m_BinomialCoefficients[n - 1][mm - 1];
      }
    }

    /// Return the binomial coefficient C_n,k. Use pre-computed values if possible, or use recursion otherwise
    int64_t BinomialCoefficient(int n, int k) const {
      if (k < 0 || k > n) return 0;
      if (k == 0 || k == n) return 1;
      if (n < m_Nmax + 1) return m_BinomialCoefficients[n][k];
      return BinomialCoefficient(n - 1, k) + BinomialCoefficient(n - 1, k - 1);
    }

    /// Returns the r-th (central) moment
    /// \param  r               Order of the moment
    /// \param  central_moment  Whether computed moment is central (true) or ordinary (false)
    double CalculateMoment(int r, bool central_moment) {
      if (central_moment)
        return GetCentralMoment(r);
      else
        return GetMoment(r);
    }

    /// Returns the sample covariance for r-th and q-th moments
    /// \param  r               Order of the first moment
    /// \param  q               Order of the sefond moment
    /// \param  central_moment  Whether the covariance is for central (true) or ordinary (false) moments
    double CalculateMomentsSampleCovariance(int r, int q, bool central_moment) {
      if (central_moment)
        return GetCentralMomentsSampleCovariance(r, q);
      else
        return GetMomentsSampleCovariance(r, q);
    }

    /// Returns the r-th order cumulant \kappa_r computed from (central) moments
    /// \param  r                     Order of the cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateCumulantFromMoments(int r, bool use_central_moments = true) {
      if (r == 1)
        return GetMean();

      auto cumulant_vs_centralmoments = JointCumulantToCentralMoments(std::vector<int>(r, 0), 1, !use_central_moments);
      double ret = 0.;
      for (const auto &term : cumulant_vs_centralmoments) {
        double tmp = static_cast<double>(term.second);
        for (const auto &multiplier : term.first) {
          tmp *= CalculateMoment(multiplier[0], use_central_moments);
        }
        ret += tmp;
      }
      return ret;
    }

    /// Returns the error estimate for the r-th order cumulant \kappa_r computed from (central) moments
    /// \param  r                     Order of the cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateCumulantErrorFromMoments(int r, bool use_central_moments = true) {
      if (r == 1)
        return GetMeanError();

      auto cumulant_vs_centralmoments = JointCumulantToCentralMoments(std::vector<int>(r, 0), 1, !use_central_moments);

      double ret = 0.;

      for (const auto &term1 : cumulant_vs_centralmoments) {
        double tmp1 = static_cast<double>(term1.second);
        for (const auto &multiplier1 : term1.first) {
          tmp1 *= CalculateMoment(multiplier1[0], use_central_moments);
        }

        for (const auto &multiplier1 : term1.first) {
          double terr1 = tmp1 / CalculateMoment(multiplier1[0], use_central_moments);

          for (const auto &term2 : cumulant_vs_centralmoments) {
            double tmp2 = static_cast<double>(term2.second);
            for (const auto &multiplier2 : term2.first) {
              tmp2 *= CalculateMoment(multiplier2[0], use_central_moments);
            }

            for (const auto &multiplier2 : term2.first) {
              double terr2 = tmp2 / CalculateMoment(multiplier2[0], use_central_moments);

              ret += terr1 * terr2 *
                     CalculateMomentsSampleCovariance(multiplier1[0], multiplier2[0], use_central_moments);
            }
          }

        }
      }

      return std::sqrt(ret);
    }

    /// Returns the cumulant ratio \kappa_r/\kappa_q computed from (central) moments
    /// \param  r                     Order of the first cumulant
    /// \param  q                     Order of the second cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateCumulantRatioFromMoments(int r, int q, bool use_central_moments = true) {
      return CalculateCumulantFromMoments(r, use_central_moments) /
             CalculateCumulantFromMoments(q, use_central_moments);
    }

    /// Returns the error estimate for the cumulant ratio \kappa_r/\kappa_q computed from (central) moments
    /// \param  r                     Order of the first cumulant
    /// \param  q                     Order of the second cumulant
    /// \param  use_central_moments   Whether to use central (true) or ordinary (false) moments for the calculation
    double CalculateCumulantRatioErrorFromMoments(int r, int q, bool use_central_moments = true) {
      if (q == 1 && use_central_moments)
        return GetScaledCumulantError(r);

      if (r == 1 && use_central_moments) {
        return GetCumulantRatio(r, q) * GetCumulantRatio(r, q) * GetScaledCumulantError(q);
      }

      double numerator = CalculateCumulantFromMoments(r, use_central_moments);
      double denominator = CalculateCumulantFromMoments(q, use_central_moments);

      auto cumulant_vs_centralmoments_numerator = JointCumulantToCentralMoments(std::vector<int>(r, 0), 1,
                                                                                !use_central_moments);
      auto cumulant_vs_centralmoments_denominator = JointCumulantToCentralMoments(std::vector<int>(q, 0), 1,
                                                                                  !use_central_moments);

      double ret = 0.0;

      // First loop over the terms in the numerator
      for (const auto &term1num : cumulant_vs_centralmoments_numerator) {
        double tmp1 = static_cast<double>(term1num.second);
        for (const auto &multiplier1num : term1num.first) {
          tmp1 *= CalculateMoment(multiplier1num[0], use_central_moments);
        }

        for (const auto &multiplier1num : term1num.first) {
          double terr1 = tmp1 / CalculateMoment(multiplier1num[0], use_central_moments);

          // Second term from the numerator
          for (const auto &term2num : cumulant_vs_centralmoments_numerator) {
            double tmp2 = static_cast<double>(term2num.second);
            for (const auto &multiplier2num : term2num.first) {
              tmp2 *= CalculateMoment(multiplier2num[0], use_central_moments);
            }

            for (const auto &multiplier2num : term2num.first) {
              double terr2 = tmp2 / CalculateMoment(multiplier2num[0], use_central_moments);

              ret += (terr1 / denominator) * (terr2 / denominator)
                     * CalculateMomentsSampleCovariance(multiplier1num[0], multiplier2num[0], use_central_moments);
            }
          }

          // Second term from the denominator, double the contribution due to symmetry
          for (const auto &term2den : cumulant_vs_centralmoments_denominator) {
            double tmp2 = static_cast<double>(term2den.second);
            for (const auto &multiplier2den : term2den.first) {
              tmp2 *= CalculateMoment(multiplier2den[0], use_central_moments);
            }

            for (const auto &multiplier2den : term2den.first) {
              double terr2 = tmp2 / CalculateMoment(multiplier2den[0], use_central_moments);

              ret += 2. * (terr1 / denominator) * (terr2 * numerator * (-1.) / denominator / denominator)
                     * CalculateMomentsSampleCovariance(multiplier1num[0], multiplier2den[0], use_central_moments);
            }
          }

        }
      }

      // Denominator-denominator terms
      for (const auto &term1den : cumulant_vs_centralmoments_denominator) {
        double tmp1 = static_cast<double>(term1den.second);
        for (const auto &multiplier1den : term1den.first) {
          tmp1 *= CalculateMoment(multiplier1den[0], use_central_moments);
        }

        for (const auto &multiplier1den : term1den.first) {
          double terr1 = tmp1 / CalculateMoment(multiplier1den[0], use_central_moments);

          // Second term from the denominator
          for (const auto &term2den : cumulant_vs_centralmoments_denominator) {
            double tmp2 = static_cast<double>(term2den.second);
            for (const auto &multiplier2den : term2den.first) {
              tmp2 *= CalculateMoment(multiplier2den[0], use_central_moments);
            }

            for (const auto &multiplier2den : term2den.first) {
              double terr2 = tmp2 / CalculateMoment(multiplier2den[0], use_central_moments);

              ret += (terr1 * numerator * (-1.) / denominator / denominator) *
                     (terr2 * numerator * (-1.) / denominator / denominator)
                     * CalculateMomentsSampleCovariance(multiplier1den[0], multiplier2den[0], use_central_moments);
            }
          }
        }
      }

      return std::sqrt(ret);
    }

    /// Calculate the values of the ordinary and central moments from the currently accumulated observations
    void CalculateMoments() {
      if (m_MomentsComputed)
        return;

      if (m_NumberOfObservations > 0) {
        for (int i = 0; i < m_Nmax + 1; ++i)
          m_Moments[i] = m_MomentSums[i] / m_NumberOfObservations;
      } else {
        m_MomentsComputed = true;
        m_nE = 0.0;
        return;
      }

      m_CentralMoments[0] = 1.;
      m_CentralMoments[1] = 0.;
      for (int mm = 2; mm < m_Nmax + 1; ++mm) {
        m_CentralMoments[mm] = 0.;
        double tmpn = 1.;
        for (int i = 0; i <= mm; ++i) {
          m_CentralMoments[mm] += BinomialCoefficient(mm, i) * m_Moments[mm - i] * tmpn * ((i & 1) ? -1. : 1.);
          tmpn *= m_Moments[1];
        }
      }

      // if m_MomentSums are mean-shifted, need to recalculate moments from the shifted moments
      if (IsMeanShifted()) {
        auto shifted_moments = m_Moments;
        m_Moments[0] = 1.;
        for (int n = 1; n <= m_Nmax; ++n) {
          m_Moments[n] = 0.;
          double Nshift = 1.;
          for (int k = 0; k <= n; ++k) {
            m_Moments[n] += BinomialCoefficient(n, k) * shifted_moments[n - k] * Nshift;
            Nshift *= m_MeanShift;
          }
        }
      }

      m_nE = static_cast<double>(m_NumberOfObservations);

      m_MomentsComputed = true;
    }


  private:
    void UpdateMomentSumsWithObservation(double N, std::vector<double> &moment_sums) {
      if (IsMeanShifted())
        N -= m_MeanShift;

      double tN = 1.;
      for (int i = 0; i < m_Nmax + 1; ++i) {
        moment_sums[i] += tN;
        tN *= N;
      }
    }

    void UpdateMomentSumsWithObservation(int N, std::vector<double> &moment_sums) {
      UpdateMomentSumsWithObservation(static_cast<double>(N), moment_sums);
    }
  };


} // namespace SampleMoments

#endif
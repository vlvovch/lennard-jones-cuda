#ifndef SPLINEFUNCTION_H
#define SPLINEFUNCTION_H
#include <vector>
#include <algorithm>
#include <cmath>

class SplineFunction
{
public:
    std::vector< std::pair<double, double> > vals;
    SplineFunction();
    SplineFunction(const std::vector<double> &x, const std::vector<double> &y)
    {
        for(int i=0;i<x.size();++i)
        {
            vals.push_back(std::make_pair(x[i], y[i]));
        }
        std::sort(vals.begin(), vals.end());
    }
    void add_val(double x, double val)
    {
        vals.push_back(std::make_pair(x, val));
        std::sort(vals.begin(), vals.end());
    }
    double f(double arg)
    {
        int ind = 0;
        std::pair<double, double> op = std::make_pair(arg, 0.);
        std::vector< std::pair<double, double> >::iterator it = std::lower_bound(vals.begin(), vals.end(), op);
        ind = std::distance(vals.begin(), it);
        /*while (ind<vals.size() && arg>=vals[ind].first)
        {
            if (fabs(vals[ind].first-arg)<1e-10) return vals[ind].second;
            ind++;
        }*/
        if (ind==0) return vals[0].second +
            (arg - vals[0].first) *
            (vals[1].second - vals[0].second) / (vals[1].first - vals[0].first);
        if (ind==vals.size()) return vals[ind-2].second +
            (arg - vals[ind-2].first) *
            (vals[ind-1].second - vals[ind-2].second) / (vals[ind-1].first - vals[ind-2].first);
        return vals[ind-1].second +
            (arg - vals[ind-1].first) *
            (vals[ind].second - vals[ind-1].second) / (vals[ind].first - vals[ind-1].first);
    }
    double fsquare(double arg)
    {
        double ret = f(arg);
        return ret*ret;
    }
    void clear()
    {
        vals.resize(2);
        vals[0].first = 0.;
        vals[0].second = 0.;
        vals[1].first = 1.;
        vals[1].second = 0.;
    }
    void clearall()
    {
        vals.resize(0);
    }
    void fill(const std::vector<double> &x, const std::vector<double> &y)
    {
        vals.resize(0);
        for(int i=0;i<x.size();++i)
        {
            vals.push_back(std::make_pair(x[i], y[i]));
        }
        std::sort(vals.begin(), vals.end());
    }
    void setConstant(double val)
    {
        vals.resize(0);
        vals.push_back(std::make_pair(0., val));
        vals.push_back(std::make_pair(1., val));
    }
    void loadFromFile(const char *file);
};

#endif // SPLINEFUNCTION_H

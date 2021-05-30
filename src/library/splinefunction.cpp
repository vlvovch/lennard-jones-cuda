#include "splinefunction.h"
#include <fstream>

using namespace std;

SplineFunction::SplineFunction()
{
    vals.resize(0);
}

void SplineFunction::loadFromFile(const char *file)
{
    fstream fin(file, fstream::in);
    vals.resize(0);
    double tx, ty;
    while (fin >> tx >> ty)
    {
        vals.push_back(make_pair(tx,ty));
    }
    sort(vals.begin(), vals.end());
}

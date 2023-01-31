#include <bits/stdc++.h>
using namespace std;
using prob_t = long long;
using ans_t = prob_t;
const prob_t PLIST_SIZE = 2000000;

struct point
{
    prob_t x, y;
};

// Return square of distance
double distance_sq(point c1, point c2)
{
    double dx = fabs(c1.x - c2.x);
    double dy = fabs(c1.y - c2.y);
    return dx * dx + dy * dy;
}

double closest(vector<point> points){
    sort(points.begin(), points.end(), [](const auto& a, const auto& b) { return a.x < b.x; });
    double mindistsq = distance_sq(points[0], points.back());
    for (size_t i = 0; i < PLIST_SIZE; i++)
    {
        for (size_t j = i + 1; j < PLIST_SIZE; j++)
        {
            if (prob_t xdist = points[j].x - points[i].x; (xdist * xdist) > mindistsq)
                break;
            mindistsq = min(mindistsq, distance_sq(points[i], points[j]));
        }
    }
    return mindistsq;
}

int main()
{
    vector<point> points;
    for (prob_t S = 290797, k = 0; k < PLIST_SIZE; k++)
    {
        prob_t x = S;
        prob_t y = S = (S * S) % 50515093;
        points.push_back({x, y});
        S = (S * S) % 50515093;
    }
    double mindistsq = closest(points);
    printf("\n%.9f\n", sqrt(mindistsq));
}

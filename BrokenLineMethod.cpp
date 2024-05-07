#include <vector>
#include <list>
#include <iostream>
// для использования бесконечности
#include <limits>
#include "BrokenLineMethod.h"
// нормирование координат?
#define INFINITY std::numeric_limits<double>::infinity()

//============================================================PROTECTED============================================================

void BrokenLineMethod::calculateR(const shared_ptr<Segment>& s) {
	for (int i = 1; i < constraintsCount + 1; i++) {
		if (s->c[i] > 0) {
			s->R = INFINITY;
			return;
		}
	}
	pair<double, double> tmp = pair<double, double>(s->p1->x, s->p2->x);
	for (int i = 1; i < constraintsCount + 1; i++) {
		double x1 = s->w[i] + s->c[i] / s->L[i];
		double x2 = s->w[i] - s->c[i] / s->L[i];

		tmp = intersectSegments(tmp, pair<double, double>(x1, x2));
	}

	double x = 0;
	if (s->w[0] >= tmp.first && s->w[0] <= tmp.second) {
		x = s->w[0];
	}
	else {
		if (s->w[0] < tmp.first)
			x = s->p1->x;
		if (s->w[0] > tmp.second)
			x = s->p2->x;
	}
	s->argmin = x;
	double f_ = s->c[0] + s->L[0] * abs(x - s->w[0]);
	s->R = f_;
}

void BrokenLineMethod::calculateLowerLimit(const shared_ptr<Segment>& s) {
	for (int i = 0; i < constraintsCount + 1; i++) {
		s->c[i] = (s->p2->f[i] + s->p1->f[i]) / 2 - s->L[i] * (s->p2->x - s->p1->x) / 2;
		s->w[i] = (s->p2->x + s->p1->x) / 2 - (s->p2->f[i] - s->p1->f[i]) / (2 * s->L[i]);
	}
}

double BrokenLineMethod::findLLoc(const shared_ptr<Segment>& s, const int i) const {
	return abs(s->p2->f[i] - s->p1->f[i]) / (s->p2->x - s->p1->x);
}



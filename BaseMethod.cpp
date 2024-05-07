#include "BaseMethod.h"
#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include <fstream>
#include <algorithm>

//TODO remove
int viewIterationsCount = 14;//Итерация, когда выводятся в файл viewData.txt
// все данные для отображения графиков функций и их нижних оценок


#define INFINITY std::numeric_limits<double>::infinity()

//============================================================PUBLIC============================================================

BaseMethod::BaseMethod(shared_ptr<TProblem> _t)
	: t(_t)
	, constraintsCount(t->getConstraintsCount())
{
	points.reserve(maxIterationCount + 2);
	degree_of_increase.resize(constraintsCount + 1);
	L_glob.resize(constraintsCount + 1);
	fileNamePoints = "resultPoints.txt";
	fileOperationCharacteristics = "resultIterationInformation.txt";
}

BaseMethod::~BaseMethod() {
	segments.clear();
	points.clear();
}

void BaseMethod::setLogNeeded(const bool _logNeeded) {
	logNeeded = _logNeeded;
}

void BaseMethod::setEps(const double _eps) {
	eps = _eps;
}

void BaseMethod::setMaxIterationsCount(const int iterationsCount)
{
	maxIterationCount = iterationsCount;
	points.reserve(maxIterationCount + 2);
}

void BaseMethod::setStopCriterion(const StopCriterion criterion) {
	stopCriterion = criterion;
}

void BaseMethod::setParams(const int N, const Gamma gamma,
	                       const double v, const double d)
{
	this->N = N;
	this->gamma1 = gamma;
	this->d = d;
	this->v = v;
}

int BaseMethod::getIterationsCount() const {
	return iterationsCount;
}

void BaseMethod::calculate()
{	
	if (logNeeded == true) {
		//вывод в файл подробной информации о каждой итерации
		ofstream streamIterationInfo(fileOperationCharacteristics);
		streamIterationInfo.close();
	}

	iterationsCount = 0;
	auto p1 = make_shared<MyPoint>(0, constraintsCount + 1, t);
	auto p2 = make_shared<MyPoint>(1, constraintsCount + 1, t);
	points.push_back(p1);
	points.push_back(p2);
	auto s = make_shared<Segment>(move(p1), move(p2), constraintsCount + 1);
	segments.push_back(s);
	calculateNewSegment(s);

	for (int i = 0; i < constraintsCount + 1; i++) {
		s->L_care[i] = max(s->L_loc[i], delta);
		L_glob[i] = s->L_loc[i];
	}

	//defineAccVal(s, s->p1);
	//defineAccVal(s, s->p2);

	recalcSegmentCharacteristics(s);

	st = s;
	double Rmin = s->R;
	if (Rmin == INFINITY) {
		st->argmin = 0.5;
	}
	iterationsCount++;
	if (logNeeded == true) {
		st->print(fileOperationCharacteristics, iterationsCount);
	}
	while (true) {
		if (needStop()) break;
		insertNewPoint(st);
		iterationsCount++;
		// TODO remove
		if (iterationsCount == ::viewIterationsCount)
			printViewData();
		if (logNeeded == true) {
			st->print(fileOperationCharacteristics, iterationsCount);
		}
	}

	if (logNeeded == true) {
		//вывод в файл всех точек измерения
		ofstream streamPoints(fileNamePoints);
		for (auto& point : points) {
			streamPoints << point->x << endl;
		}
		streamPoints.close();
	}

	cout << "xMin= " << st->argmin << endl;
	if (accVal == true) {
		cout << "yMin= " << (*t)(st->argmin, 0) << endl;
	}
	else {
		cout << "yMin= " << (*t)(st->argmin, 1) << endl;
	}
}

//============================================================PROTECTED=========================================================

void BaseMethod::defineAccVal(const shared_ptr<Segment>& s, const shared_ptr<MyPoint>& p) {
	vector<double> minInPoint;
	minInPoint.resize(constraintsCount + 1);

	bool flag = true;
	double Fmax = -INFINITY;
	// максимум значений функции ограничения в точке измерения
	for (int i = 1; i < constraintsCount + 1; i++) {
		Fmax = max(Fmax, p->f[i]);
	}
	// выполняется если хотя бы одно ограничения >0 => решаем задачу при отсутствии допустимых измерений
	if (Fmax > 0)
		flag = false;
	// появилась точка допустимых измерений
	if (flag == true) {
		accVal = true;
		minDiscrep = 0;
		currDiscrep = 0;
	}
	// понижаем уровень отсечки
	if (!accVal) {
		if (Fmax < minDiscrep) {
			minDiscrep = Fmax;
		}
		currDiscrep = minDiscrep;
	}
}


bool BaseMethod::isFullRecalcNeeded() {
	bool recalcNeeded = L_glob_changed || (needChangeGamma && segments.size() % changeGammaStep == 0);
	L_glob_changed = false;
	return recalcNeeded;
}

void BaseMethod::splitSegmentWithNewPoint(const shared_ptr<Segment>& s, shared_ptr<Segment>& s1, shared_ptr<Segment>& s2) {
	auto curIt = find(segments.begin(), segments.end(), s);
	double newValue = (*curIt)->argmin;
	auto newPoint = make_shared<MyPoint>(newValue, constraintsCount + 1, t);
	
	s1 = make_shared<Segment>((*curIt)->p1, newPoint, constraintsCount + 1);
	s2 = make_shared<Segment>(newPoint, (*curIt)->p2, constraintsCount + 1);
	points.push_back(move(newPoint));

	calculateNewSegment(s1);
	defineAccVal(s1, s1->p2);
	calculateNewSegment(s2);

	auto tmpIt = segments.erase(curIt);
	segments.insert(tmpIt, { s1, s2 });
}

void BaseMethod::recalcSegmentCharacteristics(const shared_ptr<Segment>& s) {
	calculateL_new(s);
	calculateLowerLimit(s);
	//if (accVal == false) {
	//	defineAccVal(s, s->p1);
	//	defineAccVal(s, s->p2);
	//}
	calculateR(s);
}

pair<double, double> BaseMethod::intersectSegments(const pair<double, double> x, const pair<double, double> y) const {
	pair<double, double> newX = x.first <= y.first ? x : y;
	pair<double, double> newY = x.first <= y.first ? y : x;
	if (newY.first >= newX.first && newY.second <= newX.second) {
		return pair<double, double>(newY.first, newY.second);
	}

	if (newY.first > newX.second && newY.first > newX.first) {
		return pair<double, double>(-INFINITY, INFINITY);
	}

	if (newY.first >= newX.first && newY.second > newX.second) {
		return pair<double, double>(newY.first, newX.second);
	}
	return pair<double, double>(-INFINITY, INFINITY);
}


void BaseMethod::calculateL_new(const shared_ptr<Segment>& s) {
	Gamma g = gamma1;
	if (needChangeGamma && segments.size() % changeGammaStep == 0) {
		g = gamma2;
	}

	double L_over_loc = delta;
	for (size_t i = 0; i < constraintsCount + 1; i++) {
		//расчет завышенной глобальной оценки
		calculateL_glob(s, i, g.g, g.g_constraints);

		L_over_loc = delta;
		if (s->L_care[i] > delta) {
			L_over_loc = (i == 0 ? g.g_loc : g.g_loc_constraints) * s->L_care[i]; // изменяем целевую или ограничения?
		}

		double deltaX = s->p2->x - s->p1->x;
		if (deltaX <= d) { // !!! d заранее известно или берется из выражения B*(b - a)
			s->L[i] = s->L_glob[i] * (deltaX / d) + (1 - deltaX / d) * L_over_loc;
		}
		else {
			s->L[i] = s->L_glob[i];
		}
	}
}

void BaseMethod::calculateL_care(const shared_ptr<Segment>& s) {
	auto it = find(segments.begin(), segments.end(), s);
	//для расчета осторожной оценки
	for (size_t i = 0; i < constraintsCount + 1; i++) {
		auto current = it;
		auto next = current + 1;
		if (current == segments.begin()) {
			s->L_care[i] = max({ s->L_loc[i], (*next)->L_loc[i], delta });
		}
		else {
			auto prevSegment = prev(current);
			double lprev = (*prevSegment)->L_loc[i];

			if (next == segments.end()) {
				s->L_care[i] = max({ lprev, s->L_loc[i], delta });
			}
			else {
				s->L_care[i] = max({ lprev, s->L_loc[i], (*next)->L_loc[i], delta });
			}
		}
	}
}

void BaseMethod::calculateL_glob(const shared_ptr<Segment>& s, const int i, const double g, const double g_c) {
	//расчет завышенной глобальной оценки
	if (segments.size() < N) {
		if (L_glob[i] > 0) {
			s->L_glob[i] = (i == 0 ? g : g_c) * L_glob[i]; // Ограничение целевой функции или функции ограничения?
		}
		if (L_glob[i] == 0) {
			s->L_glob[i] = 1;
		}
	}
	if (segments.size() == N) {
		if (L_glob[i] > 0) {
			degree_of_increase[i] = ((i == 0 ? g : g_c) - 1) * L_glob[i]; // вычисление delta, при k = k*
		}
		else {
			degree_of_increase[i] = 1;
		}
	}
	if (segments.size() >= N) {
		if (L_glob[i] > 0) {
			s->L_glob[i] = L_glob[i] + degree_of_increase[i];
		}
		if (L_glob[i] == 0) {
			s->L_glob[i] = 1;
		}
	}
}

// Вычисление f*
void BaseMethod::calculateFminGlob(shared_ptr<MyPoint> newPoint) {
	newfMinGlob = false;
	if (newPoint->inArea && newPoint->f[0] < fMinGlob) {
		fMinGlob = newPoint->f[0];
		newfMinGlob = true;
	}
}
//============================================================PRIVATE===========================================================

bool BaseMethod::needStop() const {
	double diff1 = abs(st->p2->x - st->argmin);
	double diff2 = abs(st->p1->x - st->argmin);
	int globMinCount = t->getGlobalSolutionCount();
	double minDif = INFINITY;
	switch (stopCriterion)
	{
	case StopCriterion::CutLength:
		return/* min(abs(st->p2->x - st->argmin), abs(st->p1->x - st->argmin)) < eps
			||*/ iterationsCount >= maxIterationCount
			|| abs(st->p2->x - st->p1->x) < eps;
		break;
	case StopCriterion::GeneratorClose:
		for (int i = 0; i < globMinCount; ++i) {
			double diff = abs((*t)[i].x - st->argmin);
			minDif = diff < minDif ? diff : minDif;
		}
		return ((minDif < eps) || iterationsCount >= maxIterationCount);
		break;
	default:
		for (int i = 0; i < globMinCount; ++i) {
			double diff = abs((*t)[i].x - st->argmin);
			minDif = diff < minDif ? diff : minDif;
		}
		return ((minDif < eps)
			|| iterationsCount >= maxIterationCount);
		break;
	}
}

void BaseMethod::insertNewPoint(const shared_ptr<Segment>& s) {
	shared_ptr<Segment> s1 = nullptr;
	shared_ptr<Segment> s2 = nullptr;
	splitSegmentWithNewPoint(s, s1, s2);

	recalcAdjacentSegments(s1);

	if (isFullRecalcNeeded()) {
		for (auto segment : segments) {
			recalcSegmentCharacteristics(segment);
		}
	}
	defineNextIteration();
	//BAG
	/*for (auto& segment : segments)
	{
		cout << segment->R << endl;
	}
	cout << "@" << endl;*/
}

void BaseMethod::defineNextIteration() {
	bool isRChanged = false;
	double Rmin = INFINITY;
	for (auto& segment : segments)
	{
		if (segment->R < Rmin)
		{
			Rmin = segment->R;
			st = segment;
			isRChanged = true;
		}
	}

	if (!isRChanged) {
		double currentMaxLength = -INFINITY;
		for (auto& segment : segments)
		{
			if (segment->p2->x - segment->p1->x > currentMaxLength)
			{
				Rmin = segment->R;
				st = segment;
				currentMaxLength = segment->p2->x - segment->p1->x;
			}
		}
		st->argmin = (st->p2->x + st->p1->x) / 2; // делим по центру
	}
	else { // делим по хитрому, формула (7) из статьи
		if (abs(st->p2->x - st->argmin) < v * (st->p2->x - st->p1->x)) {
			st->argmin = st->p2->x - v * (st->p2->x - st->p1->x);
		}
		if (abs(st->p1->x - st->argmin) < v * (st->p2->x - st->p1->x)) {
			st->argmin = st->p1->x + v * (st->p2->x - st->p1->x);
		}
	}
}

void BaseMethod::recalcAdjacentSegments(const shared_ptr<Segment>& s) {
	auto s1It = find(segments.begin(), segments.end(), s);
	//### нужен пересчет характеристик для s1, s2 и соседних
	int count = 4; // сколько сегментов пересчитываем
	std::vector<shared_ptr<Segment>>::iterator tmpIt;
	if (s1It == segments.begin()) {
		tmpIt = s1It;
		count = 3;
		std::advance(s1It, 2);
		if (s1It == segments.end()) {
			count = 2;
		}
	}
	else {
		tmpIt = prev(s1It);
		std::advance(s1It, 2);
		if (s1It == segments.end()) {
			count = 3;
		}
	}

	for (int j = 0; j < count; j++) {
		calculateL_care(*tmpIt);
		recalcSegmentCharacteristics(*tmpIt);
		tmpIt++;
	}
}

void BaseMethod::calculateNewSegment(const shared_ptr<Segment>& s)
{
	for (int i = 0; i < constraintsCount + 1; i++) {
		s->L_loc[i] = findLLoc(s, i);
		if (L_glob[i] < s->L_loc[i]) 
		{
			L_glob[i] = s->L_loc[i];
			L_glob_changed = true;
		}
	}
}

double BaseMethod::findLLocMax(const int i) const {
	double cur_max = -INFINITY;
	for (auto& segment : segments)
	{
		cur_max = max(segment->L_loc[i], cur_max);
	}
	return cur_max;
}
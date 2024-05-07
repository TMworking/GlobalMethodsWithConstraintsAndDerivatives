#include <vector>
#include <list>
#include <iostream>
#include <limits>
#include "TransformMethod.h"

#define DELTA 0.05

//============================================================PROTECTED============================================================

// вычисление характеристики интервала для количества точек > 1
void TransformMethod::calculateR(const shared_ptr<Segment>& s) {
	const int pointCoint = 200; // число точек "сеточкой" на интервале

	vector<double> maxInPoint; // элемент вектора - максимум из вычисленных значений нижних оценок в 1 точке
	maxInPoint.resize(pointCoint);

	vector<double> result; // промежуточный вектор, для вычислений в одной точке
	result.resize(constraintsCount + 1);

	double alpha = s->length / (pointCoint + 1);
	double point;
	for (int i = 0; i < pointCoint; i++) {
		point = s->p1->x + alpha * (i + 1);

		double maximum = -INFINITY;
		//по нижним оценкам ограничений
		for (int p = 0; p < constraintsCount + 1; ++p) {
			if (point < s->xt1[p]) {
				//если это парабола ветвями вниз, слева
				result[p] = s->ca[p] - (s->L[p] / 2) * (point - s->wa[p]) * (point - s->wa[p]);
			}
			else if (point > s->xt2[p]) {
				//если это парабола ветвями вниз, справа
				result[p] = s->cb[p] - (s->L[p] / 2) * (point - s->wb[p]) * (point - s->wb[p]);
			}
			else if (point <= s->xt2[p] && point >= s->xt1[p]) {
				//если это парабола ветвями вверх, в центре
				result[p] = s->c[p] + (s->L[p] / 2) * (point - s->w[p]) * (point - s->w[p]);
			}
			if (p != 0) {
				maximum = max(maximum, result[p]);
			}
			else {
				maximum = max(maximum, result[0] - fMinGlob);
			}
		}
		maxInPoint[i] = maximum;
	}

	int argminpoint = 0;
	//вектор максимумов maxInPoint по всем точкам заполнен
	//минимум среди них
	double min = maxInPoint[0];
	for (int i = 1; i < pointCoint; i++) {
		if (min > maxInPoint[i]) {
			min = maxInPoint[i];
			argminpoint = i;
		}
	}
	//характеристика интервала - найденный минимум
	s->R = min;
	s->argmin = s->p1->x + alpha * (argminpoint + 1);
}

//============================================================PRIVATE============================================================
// TODO remove
void TransformMethod::printViewData()
{
	std::ofstream out("viewData1.txt", ios_base::trunc);
	int pointsCount = 100;
	//fMinGlob = 0.97881;
	//calculateFminGlob();
	std::cout << "\nfMinGlob = " << fMinGlob << std::endl;
	vector<shared_ptr<Segment>>::iterator IT = segments.begin();
	//calculateL_new((*IT)); //Вычисление смешанных локализованных оценок 
	//calculateLowerLimit((*IT));//Вычисление параметров для склеенных парабол
							   //для всех функций
	//###
	//calculateR((*IT));

	double h = 1.0 / (pointsCount-1u),x;
	for(int i=0; i<pointsCount; ++i)
	{ 
		x = i * h;
	    if (x > (*IT)->p2->x)
	   {
			++IT;
			//calculateL_new((*IT)); //Вычисление смешанных локализованных оценок 
			//calculateLowerLimit((*IT));//Вычисление параметров для склеенных парабол
							   //для всех функций
			//###
			//calculateR((*IT));
	   }
	//Вычисление верхней огибающей функций (целевая - с трансформацией):
	   TProblem& P = (*t);
		double maxFunc = P(x, 0)-fMinGlob;
		//double y = P(x, 1);
		for (int c = 1; c <= constraintsCount; ++c)
			maxFunc = std::max(maxFunc, P(x, c)); //Верхняя огибающая функций (с трансформацией)
	//Вычисление верхней огибающей минорант ограничений 
	//  и трансформированной целевой функции: 
		double maxFuncMinor = calculateMinorant((*IT), x);
		//if (x > 0.25 && x < 0.6)
			out << x << " " << maxFunc << " " << maxFuncMinor << std::endl;

	}
	out.close();
}

// TODO remove
double TransformMethod::calculateMinorant(const shared_ptr<Segment>& s, double x)
{   double result = -INFINITY; // вычисление нижней оценки одной из функций
	double maximum = -INFINITY; //значение верхней огибающей всех оценок в одной точке
		//По нижним оценкам ограничений и целевой функции
		for (int p = 0; p < constraintsCount + 1; ++p) {
			double tmp1 = s->ca[p] - (s->L[p] / 2) * (s->xt1[p] - s->wa[p]) * (s->xt1[p] - s->wa[p]);
			double tmp2 = s->c[p] + (s->L[p] / 2) * (s->xt1[p] - s->w[p]) * (s->xt1[p]  - s->w[p]);

			if (x < s->xt1[p]) {
				//если это парабола ветвями вниз, слева
				result = s->ca[p] - (s->L[p] / 2) * (x - s->wa[p]) * (x - s->wa[p]);
			}
			else if (x > s->xt2[p]) {
				//если это парабола ветвями вниз, справа
				result = s->cb[p] - (s->L[p] / 2) * (x - s->wb[p]) * (x - s->wb[p]);
			}
			else if (x <= s->xt2[p] && x >= s->xt1[p]) {
				//если это парабола ветвями вверх, в центре
				result = s->c[p] + (s->L[p] / 2) * (x - s->w[p]) * (x - s->w[p]);
			}
			if (p != 0) {
				maximum = max(maximum, result);
			}
			else {
				maximum = max(maximum, result - fMinGlob);
			}
		}
		return maximum;
}



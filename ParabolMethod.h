#pragma once
#include "BaseMethod.h"
#include <fstream>
#include <iostream>
using namespace std;
class ParabolMethod : public BaseMethod {
public:
	//конструктор
	ParabolMethod(shared_ptr<TProblem> _t, int _mode = 1) : BaseMethod(_t) {
		mode = _mode;
		fileNamePoints = "resultPointsParabol.txt";
		fileOperationCharacteristics = "resultIterationInformationParabol.txt";
		generateFibonacci(50);
	};
	virtual ~ParabolMethod() {};
protected:
	int mode; // режим работы (численно-аналитический или по сеточке)

	vector<double> Fibonacci;

	void setMode(int _mode) {
		mode = _mode;
	}

	void Test(const shared_ptr<Segment>& s);

	//вычислить характеристику интервала
	virtual void calculateR(const shared_ptr<Segment>& s) override;

	// численно-аналитический метод
	void calculateRAnalytics(const shared_ptr<Segment>& s);

	// метод перебором по сеточке
	void calculateRSections(const shared_ptr<Segment>& s);

	//вычислить параметры нижнего ограничения
	void calculateLowerLimit(const shared_ptr<Segment>& s) override;
	//вычислить m1, m2, m3, m4 - найти максимальное
	double findLLoc(const shared_ptr<Segment>& s, const int i) const override;
private:
	//Создание массива чисел Фибоначчи
	void generateFibonacci(int n = 100);

	// верхняя огибающая нижних оценок ограничений
	pair<double, double> getEnvelope(const shared_ptr<Segment>& s, double x1, double x2);

	//Метод Фибоначчи
	pair<double, double> fibonacciMethod(const shared_ptr<Segment>& s, pair<double, double> seg);

	//найти корни квадратного уравнения
	vector<double> findRoots(const double a, const double b, const double c) const;

	//найти оценку допустимого множества, где нужно искать решения
	vector<pair<double, double>> getNegAreaOfConstraints(const shared_ptr<Segment>& s) const;

	//найти оценку допустимого множества в случае отсутствия допустимых точек
	vector<pair<double, double>> getNegAreaOfConstraints(const shared_ptr<Segment>& s, double param) const;

	//найти интервал на котором парабола отрицательна (для параболы - ветви вниз)
	vector<pair<double, double>> validSetOfNegPar(const double a, const double b, const vector<double> roots) const;

	//найти интервал на котором парабола отрицательна (для параболы - ветви вверх)
	vector<pair<double, double>> validSetOfPosPar(const double a, const double b, const vector<double> roots) const;

	//найти объединение пересечений множеств v1, v2, v3
	vector<pair<double, double>> combine(const vector<pair<double, double>> v1, const vector<pair<double, double>> v2, const vector<pair<double, double>> v3) const;

	//найти минимальное значение положительной параболы на интервале a, b
	pair<double, double> minPos(const double x, const double y, const double a,
		                        const double b, const double c, const double arg_min) const;

	//найти минимальное значение отрицательной параболы на интервале a, b
	pair<double, double> minNeg(const double x, const double y, const double a, const double b, const double c) const;

	//пересечение множеств интервалов 
	vector<pair<double, double>> intersectUnions(const vector<pair<double, double>> x, const vector<pair<double, double>> y) const;

	//найти минимальное значение ограничения целевой функции на интервале p
	pair<double, double> findMin(const shared_ptr<Segment> s, const pair<double, double> p) const;

	// TODO remove
	virtual void printViewData() override;
	//Вычисление верхней огибающей минорант ограничений  
	double calculateMinorant(const shared_ptr<Segment>& s, double x);
};


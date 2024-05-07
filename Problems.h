// +++ ������ ���������� �������� ����� �������� ���������� �����������: v2 �� 12.07.2019
#pragma once
#include <vector>
#include <iostream>
#include <corecrt_math_defines.h>
#include <algorithm>
using namespace std;

enum class TypeProblem {AnyTask, TestClass};
enum class TestFuncType {Pinter,Sheckel,Gibson,Polinom};

using TGlobalSolution = struct Point { double x; double y; };

class TProblem;

class TestFunc
{
 protected:
	double Fa, Fb;//�������� �������� ����������� ������ �������� �������
	double FaNew, FbNew;//��������, � �������� ������� ������ ���� ���������
	int pointsCount{ 1000 };//����� ����� � ����� ��� ������ ����������� Min � Max
	double epsX{ 1.0e-08 }; //�������� ����������� ���������� ��������� �� �
	double epsY{ 1.0e-05 }; //�������� ��������� �������� � ����.���������
	vector<TGlobalSolution> globalSolutions; //����� ���������� ��������� �������� ������
	double minFunc, maxFunc;//������� ��������� �������
	//��������� ��������� ������ ������ ��������� ��������� �������� �������:
	static vector<TGlobalSolution> localTmpSolutions;
	static vector<unsigned int> Fib;//����� ���������
 public:
	TestFunc(double aNew, double bNew, double a, double b) 
		: Fa(a), Fb(b), FaNew(aNew), FbNew(bNew)
	{
		generateFib();
		/*if (!Fib.size())
		{ Fib.push_back(1); Fib.push_back(1);
		for (int i = 2; i < 46; ++i)
			Fib.push_back(Fib[i - 2] + Fib[i - 1]);
		}*/
	}
 virtual void paramGeneration() = 0; //��������� ���������� �������
         //���������� ������� � "������" ����������� �� [Fa, Fb]:
 virtual double calcTestFunc(double x) = 0;
    //�������� ���������� ������� � ����������� �� [FaNew, FbNew]:
	double operator()(double x)
	{	double xx = Fa + (Fb - Fa)*(x - FaNew) / (FbNew - FaNew);
		return calcTestFunc(xx);
	}
	    //������-� �����-�� ����-� � "������" ����������� �� [Fa, Fb]:
virtual double calcDerivativeTestFunc(double x) = 0;
    //�������� ���-�� ������-� ������� � �����-�� �� [FaNew, FbNew]:
	double operator()(double x, char t)
	{ if (t != 'd') return 0.0;
	   double K = (Fb - Fa) / (FbNew - FaNew);
	   double xx = Fa + (Fb - Fa)*(x - FaNew) / (FbNew - FaNew);
	   return calcDerivativeTestFunc(xx)*K;
	}
	double getBaseLeft() { return Fa; }
	double getBaseRight() { return Fb; }
	 //����������� ������������� � ������������ ��������, 
	// ������������ ����� ����� ��������� ���������� ���������:
    void scanTestFunction();
	double getMin() { return minFunc; }
	double getMax() { return maxFunc; }
	int globalSolutionsCount() { return (int) globalSolutions.size(); }
	TGlobalSolution operator[](int i)
	{	i = i < 0 ? 0 : i;
		i = (int)((size_t)i >= globalSolutions.size() ? globalSolutions.size() - 1 : i);
		return globalSolutions[i];
	}
//����� �������� �� [x0, x1] �� [0, 1] ������� ��������� � ���������� ��������� epsXabs
// � ������� (*obj)(x), ��� obj != NULL ��� (*obj2)(x) ��� obj2 != NULL;
// ��������� ���������� �� ������ � ��������� p:
static	void getLocalMin(double x0, double x1, double epsXabs,
		             TestFunc* obj, TProblem* obj2, TGlobalSolution & p);
static	double getRandomIn(double a, double b)
	{	double per = (double)rand() / (double)RAND_MAX;
		double r = a + per * (b - a);
		return r;
	}
static void generateFib()
	{
		if (!Fib.size())
		{
			Fib.push_back(1); Fib.push_back(1);
			for (int i = 2; i < 46; ++i)
				Fib.push_back(Fib[i - 2] + Fib[i - 1]);
		}
	}
 };

class PinterTask : public TestFunc {
	double xGlobMin; //����������� ���������� �������
public:
	static	double FA, FB;//�������� ������� ��������� ��� ������� �������
public:
	PinterTask(double aNew, double bNew, double a = PinterTask::FA, double b = PinterTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "������" ������� ��������� ���������
		paramGeneration();
		scanTestFunction();
	}
virtual	void paramGeneration()
	{ xGlobMin = getRandomIn(Fa, Fb);
	 //@ cout << "xGlobMin = " << xGlobMin << " - real coordinate" << endl;
	}
	double calcTestFunc(double x)
	{	double PinterFunc;
		PinterFunc = 0.0025*(x - xGlobMin)*(x - xGlobMin) 
			+ sin(x - xGlobMin)*sin(x - xGlobMin) 
			+ sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			 * sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin));
	    return PinterFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaPinterFunc;
		derivaPinterFunc = 0.005*(x - xGlobMin)
			+ 2.0*sin(x - xGlobMin)*cos(x - xGlobMin)
			+ 2.0*sin((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			* cos((x - xGlobMin) + (x - xGlobMin)*(x - xGlobMin))
			*(1.0+2.0*(x - xGlobMin));
		return derivaPinterFunc;
	}
};
class GibsonTask : public TestFunc {
	double parA[14], parB[14], par{ 0.0 };
public:
	static double FA, FB;//�������� ������� ��������� ��� ������� �������
public:
	GibsonTask(double aNew, double bNew, double a = GibsonTask::FA, double b = GibsonTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "������" ������� ��������� ���������
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{ for (int i = 0; i < 14; ++i)
		{ parA[i] = getRandomIn(0.0, 1.0);
		  parB[i] = getRandomIn(0.0, 1.0);
		}
	}
	double calcTestFunc(double x)
	{	double GibsonFunc = parB[0];
		for (int i = 1; i < 14; i++)
			GibsonFunc += parA[i] * sin(4 * M_PI*x*i) + parB[i] * cos(4 * M_PI*x*i);
		return GibsonFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaGibsonFunc{ 0.0 };
		for (int i = 1; i < 14; i++)
			derivaGibsonFunc += (4.0 * M_PI*i)*(parA[i] * cos(4 * M_PI*x*i)
			    - parB[i] * sin(4 * M_PI*x*i));
		return derivaGibsonFunc;
	}
};

class SheckelTask : public TestFunc {
	double parC[10], parK[10], parA[10];
public:
	static	double FA, FB;//�������� ������� ��������� ��� ������� ������
public:
	SheckelTask(double aNew, double bNew, double a = SheckelTask::FA, double b = SheckelTask::FB)
		: TestFunc(aNew, bNew, a, b)
	{  // [a, b] - "������" ������� ��������� ���������
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{
		for (int i = 0; i < 10; ++i)
		{
			parC[i] = getRandomIn(0.001, 2.0);
			parK[i] = getRandomIn(1.0, 9.0);
			parA[i] = getRandomIn(0.0, 10.0);
		}
	}
	double calcTestFunc(double x)
	{ double SheckelFunc{ 0.0 };

		for (int i = 0; i < 10; i++)
			SheckelFunc += - 1.0 / (parC[i] + parK[i]*(x - parA[i])*(x - parA[i]));
		return SheckelFunc;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaSheckelFunc{ 0.0 };
	  for (int i = 1; i < 10; i++)
	  { double tmp = parC[i] + parK[i] * (x - parA[i])*(x - parA[i]);
	    derivaSheckelFunc += 2.0*parK[i] * (x - parA[i]) / (tmp*tmp);
	  }
		return derivaSheckelFunc;
	}
};
class PolinomTask : public TestFunc {
	double parR[10], par{ 1.0 }, sign{ +1.0 };
	int rootCount;
public:
	static	double FA, FB;//�������� ������� ��������� ��� ������� "�������"
public:
	PolinomTask(double aNew, double bNew, int _rootCount, 
		        double a = PolinomTask::FA, double b = PolinomTask::FB)
		: TestFunc(aNew, bNew, a, b),rootCount(_rootCount)
	{  // [a, b] - "������" ������� ��������� ���������
		paramGeneration();
		scanTestFunction();
	}
	virtual	void paramGeneration()
	{ for (int i = 1; i < rootCount - 1; ++i)
		{ parR[i] = getRandomIn(Fa + 0.05*(Fb - Fa), Fb - 0.05*(Fb - Fa));
	//@	  cout << "root = " << parR[i] << "; ";
		}
    //@	cout << endl;
	  parR[0] = Fa;  parR[rootCount - 1] = Fb;
	  sign = getRandomIn(-1.0, 1.0);
	  sign = sign < 0.0 ? -1.0 : 1.0; //�������� ���������
	  int tmpPar = (rand() % 3) - 1; //��������� -1 ��� 0 ��� +1
	  switch (tmpPar)
	  { case -1: par = 0.1; break;
	    case 0:  par = 1.0; break;
		case +1: par = 10.0; break;
		default: par = 1.0;
	  }
	
	}
	double calcTestFunc(double x)
	{ double PolinomFunc{ 1.0 };
	  for (int i = 0; i < rootCount; i++)
			PolinomFunc = PolinomFunc *(x - parR[i]);
	  return sign * PolinomFunc * par;
	}
	double calcDerivativeTestFunc(double x)
	{ double derivaPolinomFunc{ 0.0 };
	  for (int j = 0; j < rootCount; ++j)
	  { double tmp;
	    tmp = 1.0;
		for (int i = 0; i < rootCount; ++i)
			if (i == j) continue;
			else tmp = tmp * (x - parR[i]);
		derivaPolinomFunc += tmp;
	  }
	 return sign * derivaPolinomFunc * par;
	}
};

class TProblem
{
	double Fa{ 0.0 }; double Fb{ 1.0 }; //������� ������
	TypeProblem typeProblem;//��� ������ (�� ������������ ������ ��� ��������� ������)
	int numberTask;//����� ������ � ����������� ������ ��� ����� ������ � �������� ������
	int constraintsCount;//���������� �����������
	vector<int> taskFunctions;//������ �������, �������� � ����������� ������
	TestFuncType typeMainConstraint, typeMinFunc;//��� �����.����������� � �����.�������
	vector<TestFunc*> testFunctions; //����� ������� �������� ������
	vector<double> constraintsValues;//������� ������������ ������� �������� ������
	double alpha{ 0.98 }; //������������� ������� ��������� ���� �����������
	double mu{ 0.25 }; //�������� ���� ���������� �����
	int tmpPointCount{ 500 };//������ ����� ��� ���������� ���� ����������� ���������

	int pointsCount{ 5000 };//����� ����� � ����� ��� ������ ����������� Min � Max
	double epsX{ 1.0e-08 }; //�������. �������� ����������� ���������� ��������� �� �
	double epsY{ 1.0e-05 }; //�������. �������� ��������� �������� � ����.���������
	vector<TGlobalSolution> globalSolutions; //����� ������� �������� ������
	double minFunc, maxFunc;//������� ��������� ������� ������� �� ���������� ���������
		//��������� ��������� ������ ������ �������� ��������� ��������� �������� ������
		// � ������������� ����������� �� [0, 1]:
	static vector<TGlobalSolution> localTmpSolutions;
	  //������� �������� ��������� ������� (������� ������ ����������)
	  // � ������������� ����������� �� [0, 1]:
	TGlobalSolution anyLocalSolution;
public:
	TProblem(TypeProblem type = TypeProblem::AnyTask, int numTask = 1)
	{
		typeProblem = type;
		numberTask = numTask;
		if (type == TypeProblem::AnyTask)
		{  initializeTask();//��������� ���� �������� ������
		}
		else
		{ // constraintsCount = constrCount;
		 //  generateTestTask();//��������� �������� ������
		}
	}
	~TProblem()
	{	for (auto obj : testFunctions) delete obj;
	    testFunctions.clear();
		taskFunctions.clear();
		globalSolutions.clear();
		constraintsValues.clear();

	}
	TProblem(const TProblem &) = delete;
	TProblem(const TProblem &&) = delete;
	TProblem & operator=(const TProblem&) = delete;
	TProblem & operator=(const TProblem&&) = delete;

	//��������� ���������� ���������
	void setParamsTestTask(TestFuncType _typeMinFunc,
		TestFuncType _typeMainConstraint, int constrCount = 2, double _mu = 0.25)
	{	typeMinFunc = _typeMinFunc; 
		typeMainConstraint = _typeMainConstraint;
		constraintsCount = constrCount;
		alpha = 0.98;
		mu = _mu;
	}
	//�������� ���������� ������� � ���������� x �� �������������� ������� [0,1]
	// numFunc = 0 - ������� �������
	double operator()(double x, int numFunc) 
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > constraintsCount ? constraintsCount : numFunc;
		double xx = Fa + (Fb - Fa)*x;//�������� � �������� ����������
		if (typeProblem == TypeProblem::AnyTask)
		      return anyFunction(xx, taskFunctions[numFunc]);
		else return testFunction( xx, numFunc) - constraintsValues[numFunc];
	}
	//�������� ���������� ����������� ������� � ���������� x �� �������������� ������� [0,1]
    // numFunc = 0 - ������� �������
	double operator()(double x, int numFunc, char t)
	{	if (t != 'd') return 0.0;
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > constraintsCount ? constraintsCount : numFunc;
		double xx = Fa + (Fb - Fa)*x;//�������� � �������� ����������
		if (typeProblem == TypeProblem::AnyTask)
			return derivaAnyFunction(xx, taskFunctions[numFunc])*(Fb - Fa);
		else return derivaTestFunction(xx, numFunc)*(Fb - Fa);
	}
	//���������� ������� ����������� ����������� � ���������� x �� ������.������� [0, 1]
	double G(double x)
	{ double gmax = (*this)(x, 1), g;
	  for (int i = 2; i <= constraintsCount; ++i)
	   {  g = (*this)(x, i);
			if (g > gmax) gmax = g;
	   }
	  return gmax;
	}
	//�������� ��� ���������� � ������������� ����������� ��������������� ������� F(x)
	// ��� ��������� ��������� ���������� �������� (���������� �� scanTestFunction()):
	double operator()(double x)
	{ double g = G(x), F = max<double>(0.0, g);
		if (g <= 0) F += (*this)(x, 0);
		else F += max<double>((*this)(x, 0), anyLocalSolution.y);
		return F;
	}
	//������� ������� ������:
	double left() const { return Fa; }
	double right() const { return Fb; }

	// ���-� ������� �� ����������� ������ ����� 
	// � "������" ����������� ������ [Fa, Fb]:
	double anyFunction(double x, int numStFunc);
	// ���������� ������� ��������������� �������� ������ � ����������� 
	// ������ [Fa, Fb], ������ "������" ����������� �������������� �������:
	double testFunction(double x, int numFunc)
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > (int)testFunctions.size()-1 ? (int)testFunctions.size() - 1 : numFunc;
		return (*(testFunctions[numFunc]))(x);
	}
	// ���-� ����������� ������� �� ����������� ������ �����
	// � "������" ����������� ������ [Fa, Fb]:
	double derivaAnyFunction(double x, int numStFunc);
	// ���������� ����������� ������� ��������������� �������� ������ � �����������
	// ������ [Fa, Fb], ������ "������" ����������� �������������� �������:
	double derivaTestFunction(double x, int numFunc)
	{
		numFunc = numFunc < 0 ? 0 : numFunc;
		numFunc = numFunc > (int)testFunctions.size() - 1 ? 
			                  (int)testFunctions.size() - 1 : numFunc;
		return (*(testFunctions[numFunc]))(x,'d');
	}
	void initializeTask(); //������������� ���������� ��������� ������

	//������� ������� �������� ������:
	//������ ����� ������ ���� ��� isConstraints == false
	// ��� isConstraints == false � ���� Fa, Fb ������ ������ ���������� ��������
	// "������" ������ �������������� �������,
	// ��� ����������� ������� ��� ��������� ������� ����������� (��� isConstraints == true)
	// ������������ ���� ������� ������ ���������� �� ���������� ������ Fa, Fb,
	// ������ �� ����� ������������ ������
	TestFunc* generateTestFunction(TestFuncType funcType, bool isConstraints, int rootCount);

	void generateTestTask();//��������� �������� ������
	int getConstraintsCount() const //���������� ���������� �����������
	{  return constraintsCount;
	}
	//���������� �������� �������� ������������ ������� � ������� numFunc ��� ������ alpha:
	double getConstraintValue (int numFunc, double alpha) const
	{
	 return testFunctions[numFunc]->getMin()
	  + (testFunctions[numFunc]->getMax() - testFunctions[numFunc]->getMin())*alpha;
	}
	//���������� �������� Alpha � ���������� ������:
	void setNewAlphaValue(double newAlpha)
	{	alpha = newAlpha;
		for(int i = 1; i <= constraintsCount;  ++i)
			constraintsValues[i] = getConstraintValue(i, alpha);
	}
	//������������ ����������� ���� ���������� ����� �� ����� �� M �����:
	double getProportionGoodDomain(int M);
	//������ ������ ��������� Alpha ��� ���������� ���� mu ����������� ���������
	// (������ ����������� ����������� mu)
	double seachAlphaValue();
	//������ ��������� ������������� ������� ��������� ��� ������� �����������:
	double getAlpha() const { return alpha; }
	//������ �������� ������������� ���� ����������� ���������:
	double getMu() const { return mu; }
	//����������� ������������� � ������������ �������� ��������, 
	// ������������ ������ ����� ��������� �������� ���������� ��������� �������� ������:
	void scanTestTask();
	double getMin() const { return minFunc; }
	double getMax() const { return maxFunc; }
	int getGlobalSolutionCount() //���������� ����.��������� �������� ������ 
	{
		return (int)globalSolutions.size();
	}
	//��������� ������ �� i-� ���������� �������� ������:
	TGlobalSolution operator[](int i)
	{
		i = i < 0 ? 0 : i;
		i = (int)((size_t)i >= globalSolutions.size() ? globalSolutions.size() - 1 : i);
		return globalSolutions[i];
	}
	TypeProblem getTypeProblem() const { return typeProblem; }
	int getTaskNumber() const { return numberTask; }
	TestFuncType getTypeMinFunc() const { return typeMinFunc;}//��� �����.�������
	TestFuncType getTypeMainConstraint() const { return typeMainConstraint;}//��� �����.�����������
	
};


#include <iostream>
#include <math.h>
#include <list>
#include <fstream>
#include "TMatrix.h"
#include <vector>
#include <tuple>
#include <string>

using namespace std;

struct MinValue
{
    int MinIndQ, MinIndR;
    double MinQ, MinR;
};

struct Poligon
{
    // a - левая граница, b - правая граница
    // Qc - значение функции в середине отрезка
    // R - нижняя оценка
    double a, b, Qc, R; 
};

struct Data
{
    TVector<double> x, d;
};

struct Parametrs 
{
    double alpha, gamma, r, h, Mu, Eto, eps;
};

//Строит линии уровня и направления расчёта
void SolveAndLevelLines(int ItCount);
//Метод Ньютона-Рафсона
int NewtonRaphsonMethod(Parametrs& Ch, TVector<double> start, int Task_Number);
//Метод градиентного поиска
int GradientSearch(Parametrs& Ch, TVector<double> start, int Task_Number);

//Метод ломаных
void PolylineMethod(int iterations);
//Минимум функции из двух значений
double MinOfTwo(double first, double second);

//Метод деления на три
void MethodDividingByThree(int iterations);
//Следующие три отрезка
void NextPoligon(double a, double b, double L, list<Poligon>& list);
//Минимум функции Q и оценки R
MinValue MinQR(list<Poligon>& q);

//Правило Армихо
tuple<double, bool> ArmihoRule(Data& v, double Mu, double alpha, double gamma);
//Правило одномерной минимизации
tuple<double, bool> OneDimMinnRule(Data& v, double h, double Mu, double Eto);
//Правило Вулфа
tuple<double, bool> WolfeRule(Data& v, double alpha, double gamma, double r, double Mu, double Eto);

//Критерий существенности убывания функции
bool MatDecF(Data& v, double t, double Mu);
//Критерий близости к минимуму по направлению
bool ProxMinDirection2(Data& v, double t, double Mu);
//Модифицированный критерий 2
bool ProxMinDirection3(Data& v, double t, double Mu);

//Метод фибоначчи
void FibonacciMethod(double a, double b, double delta, double eps);
//Метод золотого сечения
void TheGoldenRatio(double a, double b, double delta);
//Методо дихотомии
void TheDichotomyMethod(double a, double b, double delta);

double FuncValue(double x);
double FValue(TVector<double> x);
TVector<double> GradFValue(TVector<double> x);
TMatrix<double> Hessian(TVector<double> x);

//Заполнить все критерии
void FillCriteria(Data& v, double Mu, double Eto, double RightBorder);
//Построить график
void BuildGraph(int index);

//Возращает индекс минимального значения
int MIN(TVector<double> v);
//Числа фибоначчи, начиная с 1
TVector<double> FibonacciNumbers(int N);

int main()
{
    system("chcp 1251");
    system("cls");

    int count = 0;
    Data v;
    v.x = TVector<double>(2);
    v.x[0] = 1; v.x[1] = 2;
    v.d = GradFValue(v.x) * (-1);
    
    Parametrs Ch;
    Ch.alpha = 0.01, Ch.gamma = 0.2;
    Ch.r = 2, Ch.h = pow(10, -10), Ch.eps = pow(10, -6);

    double t;
    bool ItOut;

    int q, p;
    cout << "Выберите задачу:\n1 - первая задача, 2 - вторая задача, 3 - третья задача\n";
    cin >> q;
    switch (q)
    {
    case 1:
        Ch.Mu = 1 / 16.0, Ch.Eto = 1 / 8.0;
        //Строим графики критериев
        FillCriteria(v, Ch.Mu, Ch.Eto, 0.1);
        for (int i = 1; i <= 3; i++)
            BuildGraph(i);
        break;

        //Смотрим правила поиска шагового множетеля
        tie(t, ItOut) = ArmihoRule(v, Ch.Mu, Ch.alpha, Ch.gamma);
        cout << "Правило Армихо: " << t << endl;
        //---
        tie(t, ItOut) = WolfeRule(v, Ch.alpha, Ch.gamma, Ch.r, Ch.Mu, Ch.Eto);
        cout << "Правило Вулфа: " << t << endl;
        //---
        tie(t, ItOut) = OneDimMinnRule(v, Ch.h, Ch.Mu, Ch.Eto);
        cout << "Правило одномерного поиска: " << t << endl;
        //---
        Ch.Mu = pow(10, -5); Ch.Eto = pow(10, -4);
        tie(t, ItOut) = OneDimMinnRule(v, Ch.h, Ch.Mu, Ch.Eto);
        cout << "Правило аккуратного одномерного поиска: " << t << endl;

    case 2:
        Ch.Mu = 1 / 16.0, Ch.Eto = 1 / 8.0;
        cout << "Выполнить градиентный поиск:\n";
        cout << "1 - по правилу Армихо\n";
        cout << "2 - по правилу Вулфа\n";
        cout << "3 - по правилу одномерного поиска\n";
        cout << "4 - по правилу \"аккуратного\" одномерного поиска\n";
        cin >> p;
        if (p <= 3)
            count = GradientSearch(Ch, v.x, p);
        else
        {
            Ch.Mu = pow(10, -5); Ch.Eto = pow(10, -4);
            if (p == 4)
                count = GradientSearch(Ch, v.x, 3);
        }
        SolveAndLevelLines(count);
        break;
    case 3:
        Ch.Mu = 1 / 16.0, Ch.Eto = 1 / 8.0;
        cout << "Выполнить поиск методом Ньютона-Рафсона:\n";
        cout << "1 - по правилу Армихо\n";
        cout << "2 - по правилу Вулфа\n";
        cout << "3 - по правилу одномерного поиска\n";
        cout << "4 - по правилу \"аккуратного\" одномерного поиска\n";
        cin >> p;
        if (p <= 3)
            count = NewtonRaphsonMethod(Ch, v.x, p);
        else
        {
            Ch.Mu = pow(10, -5); Ch.Eto = pow(10, -4);
            if (p == 4)
                count = NewtonRaphsonMethod(Ch, v.x, 3);
        }
        //SolveAndLevelLines(count);
        break;
    }

    return 0;
}

void SolveAndLevelLines(int ItCount)
{
    ifstream file("data.txt");
    TVector<double> x(2);
    string s;
    s = "set terminal wxt\n";
    s += "set grid\n";
    s += "set contour\n";
    s += "set nosurface\n";
    s += "set view 0,0\n";
    s += "set xrange [1:2]\n";
    s += "set yrange [1:2]\n";
    s += "set isosample 100,100\n";
    s += "set cntrparam levels discrete ";
    for (int i = 0; i < ItCount; i++)
    {
        file >> x[0] >> x[1];
        s += to_string(FValue(x)) + ", ";
    }
    file >> x[0] >> x[1];
    s += to_string(FValue(x)) + "\n";
    s += "set table\n";
    s += "set output 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\dat.txt'\n";
    //s += "splot (1-x**2+y)**2+(cos(y*(x**2))+x)**3\n";
    s += "splot (1-y**2+x)**2+(sin(y*(x**2))+x)**3\n";
    s += "unset table\n";
    s += "set palette rgbformulae 33, 13, 10\n";
    s += "plot 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\data.txt' with lines lc 'red' title 'd', \\\n";
    s += "'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\dat.txt' u 1:2:3 with lines lt 1 lw 1 palette title 'LL'\n";
    FILE* gpipe = _popen("C:\\Users\\gnuplot\\bin\\gnuplot.exe -persist", "w");
    fprintf(gpipe, s.c_str());
    _pclose(gpipe);
}

int NewtonRaphsonMethod(Parametrs& Ch, TVector<double> start, int Task_Number)
{
    Data v;
    double NormaGrad = 1, tk;
    int MaxIt = 100000, count = 0;
    bool IterOut = false;

    ofstream data("data.txt");
    data << start[0] << " " << start[1] << endl;

    v.x = start;
    v.d = Hessian(v.x).ReverseMatrix() * GradFValue(v.x) * (-1);
    while ((NormaGrad >= Ch.eps) && (count <= MaxIt))
    {
        switch (Task_Number)
        {
        case 1:
            tie(tk, IterOut) = ArmihoRule(v, Ch.Mu, Ch.alpha, Ch.gamma);
            break;
        case 2:
            tie(tk, IterOut) = WolfeRule(v, Ch.alpha, Ch.gamma, Ch.r, Ch.Mu, Ch.Eto);
            break;
        case 3:
            tie(tk, IterOut) = OneDimMinnRule(v, Ch.h, Ch.Mu, Ch.Eto);
            break;
        }
        v.x = v.x + v.d * tk;
        v.d = Hessian(v.x).ReverseMatrix() * GradFValue(v.x) * (-1);
        NormaGrad = sqrt(v.d * v.d);

        //Заполнение
        data << v.x[0] << " " << v.x[1] << endl;

        count++;
        if (IterOut == true)
        {
            cout << "Выход НЕ по точности\n";
            break;
        }
    }
    cout << "||GradQ|| = " << NormaGrad << endl;
    cout << "Iterations = " << count << endl;
    cout << "Qmin = " << FValue(v.x) << endl;
    cout << "xmin = " << v.x << endl;
    return count;
}

int GradientSearch(Parametrs& Ch, TVector<double> start, int Task_Number)
{
    Data v;
    double NormaGrad = 1, tk;
    int MaxIt = 100000, count = 0;
    bool IterOut = false;

    ofstream data("data.txt");
    data << start[0] << " " << start[1] << endl;

    v.x = start;
    v.d = GradFValue(v.x) * (-1);
    while ((NormaGrad >= Ch.eps) && (count < MaxIt))
    {
        switch (Task_Number)
        {
        case 1:
            tie(tk, IterOut) = ArmihoRule(v, Ch.Mu, Ch.alpha, Ch.gamma);
            break;
        case 2:
            tie(tk, IterOut) = WolfeRule(v, Ch.alpha, Ch.gamma, Ch.r, Ch.Mu, Ch.Eto);
            break;
        case 3:
            tie(tk, IterOut) = OneDimMinnRule(v, Ch.h, Ch.Mu, Ch.Eto);
            break;
        }
        v.x = v.x + v.d * tk;
        v.d = GradFValue(v.x) * (-1);
        NormaGrad = sqrt(v.d * v.d);

        //Заполнение
        data << v.x[0] << " " << v.x[1] << endl;

        count++;
        if (IterOut == true)
        {
            cout << "Выход НЕ по точности\n";
            break;
        }
    }
    cout << "||GradQ|| = " << NormaGrad << endl;
    cout << "Iterations = " << count << endl;
    cout << "Qmin = " << FValue(v.x) << endl;
    cout << "xmin = " << v.x << endl;
    return count;
}

//Метод ломаных
void PolylineMethod(int iterations)
{
    double eps = pow(10, -4), L = 2, a = 0, b = 12;
    double Q1 = FuncValue(a), Q2 = FuncValue(b), Qopt, Rt = 0;
    double x1, x2, x1opt = a, x2opt = b, xmid;

    vector<double> R(1);
    Qopt = MinOfTwo(Q1, Q2);

    list<double> x;
    list<double> ::iterator it;
    x.push_back(a);
    x.push_back(b);
    
    int j = 0;
    while (j < iterations)
    {
        //Ищем минимум Rt - оценочного значения
        Rt = pow(10, 10);
        R = vector<double>(x.size() - 1);
        it = x.begin();
        x1 = *it; it++;
        for (int i = 0; i < x.size()-1; i++)
        {
            x2 = *it;
            R[i] = (FuncValue(x1) + FuncValue(x2)) / 2.0 - L * (x2 - x1) / 2.0;
            if (R[i] < Rt)
            {
                Rt = R[i];
                x1opt = x1;
                x2opt = x2;
            }
            x1 = x2;
            it++;
        }
        if (Qopt - Rt < eps)
            break;
        //Нашли Xk+1
        xmid = (x1opt + x2opt) / 2 - (FuncValue(x2opt) - FuncValue(x1opt)) / (2.0 * L);
        //Добавили точку в список
        x.push_back(xmid);
        //Отсортировали список
        x.sort();
        //Нашли функию в новой точке
        Q1 = FuncValue(xmid);
        //Обновили оптимальное значение функции
        Qopt = MinOfTwo(Q1, Qopt);
        j++;
    }
    cout << "Qopt = " << Qopt << endl;
    cout << "Rt = " << Rt << endl;
}

double MinOfTwo(double first, double second)
{
    if (first < second)
        return first;
    else return second;
}

//Метод деления на три
void MethodDividingByThree(int iterations)
{
    ofstream file("dat.txt");
    double a = 0, b = 12, L = 2, eps = pow(10, -4), c;

    Poligon p;
    MinValue v;
    v.MinR = 10; v.MinQ = 20;

    list<Poligon> Beta;
    list<Poligon> ::iterator it;

    int j = 0;
    while (j < iterations)
    {
        //Нашли текущие отрезки
        NextPoligon(a, b, L, Beta);
        //Ищем оптимальные значения на шаге
        v = MinQR(Beta);
        //Проверяем критерий
        if ((v.MinQ - v.MinR) < eps)
            break;
        //Удаляем ненужный блок
        it = Beta.begin();
        advance(it, v.MinIndR);
        p = *it;
        a = p.a; b = p.b;
        Beta.erase(it);
        j++;
    }
    cout << "R = " << v.MinR << endl;
    cout << "Q = " << v.MinQ << endl;
}

void NextPoligon(double a, double b, double L, list<Poligon> &list)
{
    double c;
    vector<Poligon> P(3);
    P[0].a = a; P[0].b = a + (b - a) / 3.0;
    P[1].a = P[0].b; P[1].b = a + 2 * (b - a) / 3.0;
    P[2].a = P[1].b; P[2].b = b;

    for (int i = 0; i < 3; i++)
    {
        c = (P[i].a + P[i].b) / 2;
        P[i].Qc = FuncValue(c);
        P[i].R = P[i].Qc - L * (P[i].b - P[i].a) / 2;
        list.push_back(P[i]);
    }
}

//Минимум R и Q у полигона
MinValue MinQR(list<Poligon>& q)
{
    list<Poligon> ::iterator it = q.begin();
    Poligon p = *it; it++;

    MinValue v; 
    v.MinIndQ = 0; v.MinIndR = 0;
    v.MinQ = p.Qc; v.MinR = p.R;

    for (int i = 1; i < q.size(); i++)
    {
        p = *it;
        //Ищем минимум функции Q
        if (p.Qc < v.MinQ)
        {
            v.MinQ = p.Qc;
            v.MinIndQ = i;
        }
        //Ищем минимум оценочного значения R
        if (p.R < v.MinR)
        {
            v.MinR = p.R;
            v.MinIndR = i;
        }
        it++;
    }

    return v;
}

//---

//Правило Армихо
tuple<double, bool> ArmihoRule(Data& v, double Mu, double alpha, double gamma)
{
    double t;
    bool IterOut = false;
    t = alpha;

    int count = 0, MaxIt = 10000;
    while ((!MatDecF(v, t, Mu)) && (count < MaxIt))
    {
        t = gamma * t;
    }
    if (count >= MaxIt)
        IterOut = true;
    return make_tuple(t, IterOut);
}

//Метод одномерной оптимизации
tuple<double, bool> OneDimMinnRule(Data& v, double h, double Mu, double Eto)
{
    bool IterOut = false;
    double t = 0, FiT = FValue(v.x + v.d * t);
    double _t = h, Fi_T = FValue(v.x + v.d * _t);

    while (Fi_T <= FiT)
    {
        t = _t; FiT = Fi_T;
        h *= 2;
        _t = t + h; Fi_T = FValue(v.x + v.d * _t);
    }

    if (fabs(_t - h) <= pow(10, -14))
    {
        IterOut = true;
        return make_tuple(_t, IterOut);
    }
    //Метод золотого сечения
    int count = 0, s = 0, MaxIt = 10000;
    //0-x, 1-y
    TVector<double> Q(2);
    double x, y, tau;
    double a = 0, b = _t;
    tau = (sqrt(5) - 1) / 2.0;
    y = a + (b - a) * tau; Q[1] = FValue(v.x + v.d * y);
    x = b - (b - a) * tau; Q[0] = FValue(v.x + v.d * x);

    t = x;
    if (MatDecF(v, t, Mu) && ProxMinDirection2(v, t, Eto))
        return make_tuple(t, IterOut);
    t = y;
    while (!MatDecF(v, t, Mu) || !ProxMinDirection2(v, t, Eto))
    {
        s = MIN(Q);
        switch (s)
        {
        case 0:
            b = y; y = x; Q[1] = Q[0];
            x = b - (b - a) * tau;
            Q[0] = FValue(v.x + v.d * x);
            t = x;
            break;
        case 1:
            a = x; x = y; Q[0] = Q[1];
            y = a + (b - a) * tau;
            Q[1] = FValue(v.x + v.d * y);
            t = y;
            break;
        }
        count++;
        if (count >= MaxIt)
        {
            IterOut = true;
            break;
        }
    }
    return make_tuple(t, IterOut);
}

//Правило Вулфа
tuple<double, bool> WolfeRule(Data& v, double alpha, double gamma, double r, double Mu, double Eto)
{
    bool IterOut = false;

    double t = alpha;
    double alpha1 = 0, alpha2 = 0;

    int count = 0, MaxIt = 10000;
    while ((!MatDecF(v, t, Mu) || !ProxMinDirection3(v, t, Eto)) && (count < MaxIt))
    {
        while ((!MatDecF(v, t, Mu)) && (count < MaxIt))
        {
            alpha2 = t;
            t = alpha1 * gamma + alpha2 * (1 - gamma);
            count++;
        }
        while ((!ProxMinDirection3(v, t, Eto)) && (count < MaxIt))
        {
            alpha1 = t;
            if (alpha2 == 0)
            {
                t = alpha1 * r;
                break;
            }
            else
                t = alpha1 * gamma + alpha2 * (1 - gamma);
            count++;
        }
        count++;
    }
    if (count >= MaxIt)
        IterOut = true;
    return make_tuple(t, IterOut);
}

//---

//Критерий существенности убывания фукции
bool MatDecF(Data& v, double t, double Mu)
{
    bool FeedBack = false;
    TVector<double> Grad_x = GradFValue(v.x);
    if (t >= 0)
        if (FValue(v.x + v.d * t) <= FValue(v.x) + Mu * (Grad_x * v.d) * t)
            FeedBack = true;
    return FeedBack;
}

//Критерий близости к минимуму по направлению 2
bool ProxMinDirection2(Data& v, double t, double Eto)
{
    bool FeedBack = false;
    TVector<double> Grad_xd = GradFValue(v.x + v.d * t), Grad_x = GradFValue(v.x);
    if (t >= 0)
        if (fabs(Grad_xd * v.d) <= Eto * fabs(Grad_x * v.d))
            FeedBack = true;
    return FeedBack;
}

//Критерий близости к минимуму по направлению 3
bool ProxMinDirection3(Data& v, double t, double Eto)
{
    bool FeedBack = false;
    TVector<double> Grad_xd = GradFValue(v.x + v.d * t), Grad_x = GradFValue(v.x);
    if (t >= 0)
        if (Grad_xd * v.d >= Eto * (Grad_x * v.d))
            FeedBack = true;
    return FeedBack;
}

//---

void FibonacciMethod(double a, double b, double delta, double eps)
{
    int s = 0, N = 1, k;
    //0-x, 1-y
    TVector<double> Q(2);
    double x, y;
    //---
    double F1 = 1, F2 = 1, temp;
    while ((b - a) / F2 + eps > delta)
    {
        temp = F1;
        F1 = F2;
        F2 = temp + F2;
        N++;
    }
    TVector<double> F = FibonacciNumbers(N + 1);
    x = b - (b - a) * F[N - 1] / F[N];
    y = a + (b - a) * F[N - 1] / F[N];
    Q[0] = FValue(x); Q[1] = FValue(y);
    k = 2;
    for (int k = 2; k <= N; k++)
    {
        s = MIN(Q);
        switch (s)
        {
        case 0:
            if (k <= N - 2)
            {
                b = y; y = x; Q[1] = Q[0];
                x = y - (y - a) * F[N - k - 2] / F[N - k];
                Q[0] = FValue(x);
            }
            else
                if (k == N - 1)
                {
                    b = y; y = x; Q[1] = Q[0];
                    x = y - eps;
                    Q[0] = FValue(x);
                }
                else
                    b = y;
            break;
        case 1:
            if (k <= N - 2)
            {
                a = x; x = y; Q[0] = Q[1];
                y = x + (b - x) * F[N - k - 2] / F[N - k];
                Q[1] = FValue(y);
            }
            else
                if (k == N - 1)
                {
                    a = x; x = y; Q[0] = Q[1];
                    y = x + eps;
                    Q[1] = FValue(y);
                }
                else
                    a = x;
            break;
        }
    }
    cout << a << " " << b << endl;
    cout << N - 1 << endl;
}

void TheGoldenRatio(double a, double b, double delta)
{
    int count = 0, s = 0;
    //0-x, 1-y
    TVector<double> Q(2);
    double x, y, tau;
    tau = 0.62;
    y = a + (b - a) * tau; Q[1] = FValue(y);
    x = b - (b - a) * tau; Q[0] = FValue(x);
    while ((b - a) > delta)
    {
        s = MIN(Q);
        switch (s)
        {
        case 0:
            b = y; y = x; Q[1] = Q[0];
            x = b - (b - a) * tau;
            Q[0] = FValue(x);
            break;
        case 1:
            a = x; x = y; Q[0] = Q[1];
            y = a + (b - a) * tau;
            Q[1] = FValue(y);
            break;
        }
        count++;
    }
    cout << a << " " << b << endl;
    cout << count << endl;
}

void TheDichotomyMethod(double a, double b, double delta)
{
    int count = 0, s;
    //0-x, 1-y, 2-c
    TVector<double> Q(3);
    double c, x, y;
    c = (a + b) / 2.0;
    Q[2] = FValue(c);
    while ((b - a) > delta)
    {
        x = (a + c) / 2.0; y = (c + b) / 2.0;
        Q[0] = FValue(x); Q[1] = FValue(y);
        s = MIN(Q);
        switch (s)
        {
        case 0:
            b = c; c = x; Q[2] = Q[0];
            break;
        case 1:
            a = c; c = y; Q[2] = Q[1];
            break;
        case 2:
            a = x; b = y;
            break;
        }
        count++;
    }
    cout << a << " " << b << endl;
    cout << count << endl;
}

//---

double FuncValue(double x)
{
    return sin(x);
}

double FValue(TVector<double> x)
{
    //return pow(x[0], 2) + pow(x[1], 2) + 3;
    //return pow((1 - pow(x[0], 2) + x[1]), 2) + pow(cos(pow(x[0], 2) * x[1]) + x[0], 3);
    return pow((1 - pow(x[1], 2) + x[0]), 2) + pow(sin(pow(x[0], 2) * x[1]) + x[0], 3);
}

TVector<double> GradFValue(TVector<double> x)
{
    double temp1, temp2;
    /*temp1 = 1 - pow(x[0], 2) + x[1];
    temp2 = pow(cos(pow(x[0], 2) * x[1]) + x[0], 2);*/
    temp1 = 1 - pow(x[1], 2) + x[0];
    temp2 = pow(sin(pow(x[0], 2) * x[1]) + x[0], 2);

    TVector<double> res(2);
    /*res[0] = 2 * x[0];
    res[1] = 2 * x[1];*/
    /*res[0] = -4 * x[0] * temp1 + 3 * (1 - 2 * x[0] * x[1] * sin(pow(x[0], 2) * x[1])) * temp2;
    res[1] = 2 * temp1 - 3 * pow(x[0], 2) * sin(pow(x[0], 2) * x[1]) * temp2;*/
    res[0] = 2 * temp1 + 3 * (1 + 2 * x[0] * x[1] * cos(pow(x[0], 2) * x[1])) * temp2;
    res[1] = -4 * x[1] * temp1 + 3 * pow(x[0], 2) * cos(pow(x[0], 2) * x[1]) * temp2;
    return res;
}

TMatrix<double> Hessian(TVector<double> x)
{
    TMatrix<double> res(2);
    double tc, ts, tcx, txy;
    tc = cos(pow(x[0], 2) * x[1]);
    ts = sin(pow(x[0], 2) * x[1]);
    tcx = tc + x[0];
    txy = 1 - pow(x[0], 2) + x[1];
    /*res[0][0] = 2;
    res[1][0] = 0;
    res[0][1] = 0;
    res[1][1] = 2;*/
    res[0][0] = 3 * pow(tcx, 2) * (-4 * pow(x[0], 2) * pow(x[1], 2) * tc - 2 * x[1] * ts)
        - 4 * txy + 6 * pow((1 - 2 * x[0] * x[1] * ts), 2) * tcx + 8 * pow(x[0], 2);
    res[0][1] = -6 * pow(x[0], 2) * ts * (1 - 2 * x[0] * x[1] * ts) * tcx
        + 3 * pow(tcx, 2) * (-2 * x[0] * ts - 2 * pow(x[0], 3) * x[1] * tc) - 4 * x[0];
    res[1][0] = res[0][1];
    res[1][1] = -3 * pow(x[0], 4) * tc * pow(tcx, 2) + 6 * pow(x[0], 4) * pow(ts, 2) * tcx + 2;
    return res;
}

//---

void FillCriteria(Data& v, double Mu, double Eto, double RightBorder)
{
    ofstream Section("Section.txt"), P1("P1.txt"), P2("P2.txt"), P3("P3.txt");

    double h = pow(10, -6), t = 0;
    double temp;
    while (t <= RightBorder)
    {
        temp = FValue(v.x + v.d * t);
        Section << t << " " << temp << endl;
        if (MatDecF(v, t, Mu))
            P1 << t << " " << temp << endl;
        if (ProxMinDirection2(v, t, Eto))
            P2 << t << " " << temp << endl;
        if (ProxMinDirection3(v, t, Eto))
            P3 << t << " " << temp << endl;
        t += h;
    }
}

void BuildGraph(int index)
{
    FILE* gpipe = _popen("C:\\Users\\gnuplot\\bin\\gnuplot.exe -persist", "w");
    fprintf(gpipe, "set terminal wxt persist\n");
    fprintf(gpipe, "set grid\n");
    string s;
    switch (index)
    {
    case 1:
        s = "plot 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\P1.txt'";
        s += " with points lc 'red' title ' ', ";
        s += "'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\Section.txt'";
        s += " with lines lc 'blue' title ' '\n";
        break;
    case 2:
        s = "plot 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\P2.txt'";
        s += " with points lc 'red' title ' ', ";
        s += "'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\Section.txt'";
        s += " with lines lc 'blue' title ' '\n";
        break;
    case 3:
        s = "plot 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\P3.txt'";
        s += " with points lc 'red' title ' ', ";
        s += "'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\Section.txt'";
        s += " with lines lc 'blue' title ' '\n";
        break;
    case 4:
        s = "plot 'C:\\Users\\HOME\\source\\repos\\Method_Optimisation\\Method_Optimisation\\data.txt'";
        s += " with lines lc 'red' title ' '\n";
        break;

    }
    fprintf(gpipe, s.c_str());
    _pclose(gpipe);
}

//---

int MIN(TVector<double> v)
{
    int n = v.Size(), s = 0;
    double min = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] < min)
        {
            min = v[i];
            s = i;
        }
    return s;
}

TVector<double> FibonacciNumbers(int N)
{
    TVector<double> res(N);
    if (N == 1)
    {
        res[0] = 1;
        return res;
    }
    if (N > 1)
    {
        res[0] = 1; res[1] = 1;
        if (N > 2)
            for (int i = 2; i < N; i++)
                res[i] = res[i - 1] + res[i - 2];
    }
    return res;
}
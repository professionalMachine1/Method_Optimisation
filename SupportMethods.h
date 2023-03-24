//#pragma once
//
//
////Метод ломаных
//void PolylineMethod()
//{
//    double eps = pow(10, -4), L = 2, a = 0, b = 12;
//    double Q1 = FuncValue(a), Q2 = FuncValue(b), Qopt, Rt;
//    double x1, x2, x1opt = a, x2opt = b, xmid;
//    vector<double> R(1);
//    Qopt = MinOfTwo(Q1, Q2);
//    list<double> x;
//    list<double> ::iterator it;
//    x.push_back(a);
//    x.push_back(b);
//
//    ofstream file("dat.txt");
//    int j = 0;
//    while (j < 16)
//    {
//        //Строим текущее положение дел
//        it = x.begin();
//        x1 = *it; it++;
//        file << x1 << " " << FuncValue(x1) << endl;
//        for (int i = 1; i < x.size(); i++)
//        {
//            x2 = *it;
//            Q1 = FuncValue(x1); Q2 = FuncValue(x2);
//            xmid = (x1 + x2) / 2.0 - (Q2 - Q1) / (2.0 * L);
//            Rt = (Q1 + Q2) / 2.0 - L * (x2 - x1) / 2.0;
//            file << xmid << " " << Rt << endl;
//            file << x2 << " " << Q2 << endl;
//            x1 = x2;
//            it++;
//        }
//        file << "\n\n";
//
//        Rt = pow(10, 10);
//        R = vector<double>(x.size() - 1);
//        it = x.begin();
//        x1 = *it; it++;
//        for (int i = 0; i < x.size() - 1; i++)
//        {
//            x2 = *it;
//            R[i] = (FuncValue(x1) + FuncValue(x2)) / 2.0 - L * (x2 - x1) / 2.0;
//            if (R[i] < Rt)
//            {
//                Rt = R[i];
//                x1opt = x1;
//                x2opt = x2;
//            }
//            x1 = x2;
//            it++;
//        }
//        if (Qopt - Rt < eps)
//            break;
//        //Нашли Xk+1
//        xmid = (x1opt + x2opt) / 2 - (FuncValue(x2opt) - FuncValue(x1opt)) / (2.0 * L);
//        //Добавили точку в список
//        x.push_back(xmid);
//        //Отсортировали список
//        x.sort();
//        //Нашли функию в новой точке
//        Q1 = FuncValue(xmid);
//        //Обновили оптимальное значение функции
//        Qopt = MinOfTwo(Q1, Qopt);
//        j++;
//    }
//    cout << "Qopt = " << Qopt << endl;
//    cout << "Rt = " << Rt << endl;
//}
//
//double MinOfTwo(double first, double second)
//{
//    if (first < second)
//        return first;
//    else return second;
//}
//
////Метод деления на три
//void MethodDividingByThree()
//{
//    ofstream file("dat.txt");
//    double a = 0, b = 12, L = 2, eps = pow(10, -4), c;
//
//    Poligon p;
//    MinValue v;
//    v.MinR = 10; v.MinQ = 20;
//
//    list<Poligon> Beta;
//    list<Poligon> ::iterator it;
//
//    list<double> Points;
//    list<double> ::iterator ItP;
//    double temp;
//
//    int j = 0;
//    while (j < 16)
//    {
//        //Нашли текущие отрезки
//        NextPoligon(a, b, L, Beta);
//        //Добавляем их в файл
//        it = Beta.begin();
//        for (int i = 0; i < Beta.size(); i++)
//        {
//            p = *it;
//            Points.push_back(p.a);
//            Points.push_back(p.b);
//            it++;
//        }
//        //---
//        Points.unique();
//        Points.sort();
//        ItP = Points.begin();
//        a = *ItP; ItP++;
//        file << a << " " << FuncValue(a) << endl;
//        for (int i = 1; i < Points.size(); i++)
//        {
//            b = *ItP;
//            c = (a + b) / 2;
//            file << c << " " << FuncValue(c) - L * (b - a) / 2.0 << endl;
//            file << b << " " << FuncValue(b) << endl;
//            a = b;
//            ItP++;
//        }
//        file << "\n\n";
//        Points.clear();
//
//        //Ищем оптимальные значения на шаге
//        v = MinQR(Beta);
//        //Проверяем критерий
//        if ((v.MinQ - v.MinR) < eps)
//            break;
//        //удаляем ненужный блок
//        it = Beta.begin();
//        advance(it, v.MinIndR);
//        p = *it;
//        a = p.a; b = p.b;
//        Beta.erase(it);
//        j++;
//    }
//    cout << "R = " << v.MinR << endl;
//    cout << "Q = " << v.MinQ << endl;
//}
//
//void NextPoligon(double a, double b, double L, list<Poligon>& list)
//{
//    vector<Poligon> P(3);
//    P[0].a = a; P[0].b = a + (b - a) / 3.0; P[0].c = (P[0].a + P[0].b) / 2;
//    P[1].a = P[0].b; P[1].b = a + 2 * (b - a) / 3.0; P[1].c = (P[1].a + P[1].b) / 2;
//    P[2].a = P[1].b; P[2].b = b; P[2].c = (P[2].a + P[2].b) / 2;
//
//    for (int i = 0; i < 3; i++)
//    {
//        P[i].Qc = FuncValue(P[i].c);
//        P[i].R = P[i].Qc - L * (P[i].b - P[i].a) / 2;
//        list.push_back(P[i]);
//    }
//}
//
////Минимум R и Q у полигона
//MinValue MinQR(list<Poligon>& q)
//{
//    list<Poligon> ::iterator it = q.begin();
//    Poligon p = *it; it++;
//
//    MinValue v;
//    v.MinIndQ = 0; v.MinIndR = 0;
//    v.MinQ = p.Qc; v.MinR = p.R;
//
//    //Ищем минимум функции Q
//    for (int i = 1; i < q.size(); i++)
//    {
//        p = *it;
//        if (p.Qc < v.MinQ)
//        {
//            v.MinQ = p.Qc;
//            v.MinIndQ = i;
//        }
//        it++;
//    }
//
//    it = q.begin();
//    it++;
//    //Ищем минимум параметра R
//    for (int i = 1; i < q.size(); i++)
//    {
//        p = *it;
//        if (p.R < v.MinR)
//        {
//            v.MinR = p.R;
//            v.MinIndR = i;
//        }
//        it++;
//    }
//
//    return v;
//}
//

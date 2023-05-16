## Летняя практика. Обработка результатов эксперимента. Метод наименьших квадратов.

Золотухина П. В.
Группа 23/3
Вариант 4

### Задача
В результате эксперимента была определена некоторая табличная зависимость.
С помощью метода наименьших квадратов определить линию регрессии, рассчитать коэффициент корреляции, подобрать функциональную зависимость заданного вида, вычислить коэффициент регрессии. Построить график экспериментальной зависимости, линию регрессии и график подобранной зависимости. Определить суммарную квадратичную ошибку и среднюю ошибку для линии регрессии и подобранной функциональной зависимости. Написать программу на языке С(С++) для решения задачи. При необходимости напишите функцию решения системы линейных алгебраических уравнений.

![image](https://github.com/pococonut/summer_practice/assets/114181600/efae3748-f058-45bb-b218-f240a9ae8771)

#### Программа
```
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// Вычисление коэффициентов аппроксимирующей прямой
void getApprox(float (&x)[2][7], double *a, double *b, int n) {
    double sumx = 0;
    double sumy = 0;
    double sumx2 = 0;
    double sumxy = 0;
    for (int i = 0; i<n; i++) {
        sumx += x[0][i];
        sumy += x[1][i];
        sumx2 += x[0][i] * x[0][i];
        sumxy += x[0][i] * x[1][i];
    }
    *a = (n*sumxy - (sumx*sumy)) / (n*sumx2 - sumx*sumx);
    *b = (sumy - *a*sumx) / n;
    return;
}
int main() {

    setlocale(LC_ALL, "RU");
    double  a, b;
    const unsigned int DIM1 = 2;
    const unsigned int DIM2 = 7;

    float arr[DIM1][DIM2] = {
        { 0.2, 0.7, 1.2, 1.70, 2.2, 2.7, 3.2 },
        { 2.3198, 2.8569, 3.5999, 4.4357, 5.5781, 6.9459, 8.6621 },
    };

    cout << "x     y" << endl;
    for (int i = 0; i < DIM2; i++)
        cout << arr[0][i] << " - " << arr[1][i] << endl;
    getApprox(arr, &a, &b, DIM2);

    cout << endl << "a = " << a << " " << "b = " << b << endl;
    cout << endl << "Уравнение регрессии имеет вид: " << "y = " << b << " + " << a << "*x" << endl;

    //Формирование массивов для изображения линии регрессии на графике.
    float y2[DIM2];
    for (int i = 0; i < DIM2; i++)
        y2[i] = b + a * arr[0][i];

    // коэффцент корреляции
    float sumx = 0, sumy = 0, srx = 0, sry = 0;

    for (int i = 0; i < DIM2; i++) {
        sumx += arr[0][i];
        sumy += arr[1][i];
    }

    srx = sumx / DIM2;
    sry = sumy / DIM2;

    float SR[DIM1][DIM2] = {  //отклонения от среднего арифметического
        { .0, .0, .0, .0, .0, .0, .0 },
        { .0, .0, .0, .0, .0, .0, .0 },
    };

    float SR2[DIM1][DIM2] = {
        { .0, .0, .0, .0, .0, .0, .0 },
        { .0, .0, .0, .0, .0, .0, .0 },
    };

    float SRxy[DIM2] = { .0, .0, .0, .0, .0, .0, .0 };

    for (int i = 0; i < DIM2; i++)
    {
        SR[0][i] = arr[0][i] - srx;
        SR[1][i] = arr[1][i] - sry;
    }
    SR[0][3] = 0;
    
    for (int i = 0; i < DIM1; i++)
        for (int j = 0; j < DIM2; j++)
            SR2[i][j] = SR[i][j] * SR[i][j];

    for (int i = 0; i < DIM2; i++)
        SRxy[i] = SR[0][i] * SR[1][i];

    // сумма произведений отклонений
    float sum_pr_ot_xy = 0;
    for (int i = 0; i < DIM2; i++)
        sum_pr_ot_xy += SRxy[i];

    // суммы квадратов отклонений по x и y
    float sum_kv_ot_x = 0;
    float sum_kv_ot_y = 0;
    for (int j = 0; j < DIM2; j++) {
        sum_kv_ot_x += SR2[0][j];
        sum_kv_ot_y += SR2[1][j];
    }

    //коэфф кореляции
    float kor = 0;
    kor = sum_pr_ot_xy / (sqrt(sum_kv_ot_x * sum_kv_ot_y));

    cout << endl << "Коэффицент кореляции: " << kor << endl;
    cout << endl;

    float Y_reg[DIM2], SSE = 0, MAE = 0;

    for (int i = 0; i < DIM2; i++)
        Y_reg[i] = b + a * (arr[0][i]);
    
    for (int i = 0; i < DIM2; i++)
        SSE += pow(arr[1][i] - Y_reg[i], 2);

    cout << "Суммарная квадратичная ошибка для линии регрессии SSE: " << SSE << endl;

    MAE = SSE / DIM2;

    cout << endl << "Средняя ошибка для линии регрессии MAE: " << MAE << endl;


    //  Подбор параметров функции y = ax be cx
    // считаем суммы

    float sum_lnx = 0, sum_x = 0, sum_y = 0, sum_lnx2 = 0, sum_xlnx = 0, sum_ylnx = 0, sum_x2 = 0, sum_yx = 0;

    for (int i = 0; i < DIM2; i++)
    {
        sum_lnx += log(arr[0][i]);
        sum_x += arr[0][i];
        sum_y += arr[1][i];
        sum_lnx2 += pow(log(arr[0][i]), 2);
        sum_xlnx += arr[0][i] * log(arr[0][i]);
        sum_ylnx += arr[1][i] * log(arr[0][i]);
        sum_x2 += pow(arr[1][i], 2);
        sum_yx += arr[0][i] * arr[1][i];

    }

    // решение системы методом гаусса

    const unsigned int m = 3;

    double d, s, a1[m + 1][m + 1], x[m + 1];

    double M[m + 1][m+1] = { {0,0,0,0}, {0,7, sum_lnx, sum_x},
                           {0, sum_lnx, sum_lnx2, sum_xlnx},
                           {0, sum_x, sum_xlnx, sum_x2} };

    double B[m + 1] = { 0, sum_y, sum_ylnx, sum_yx };


    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            a1[i][j] = M[i][j];
        }
    }
    cout << endl;


    for (int k = 1; k <= m; k++) // прямой ход
    {
        for (int j = k + 1; j <= m; j++)
        {
            d = M[j][k] / M[k][k]; // формула (1)

            for (int i = k; i <= m; i++)
            {
                M[j][i] = M[j][i] - d * M[k][i]; // формула (2)
            }
            B[j] = B[j] - d * B[k]; // формула (3)
        }
    }

    for (int k = m; k >= 1; k--) // обратный ход
    {
        d = 0;
        for (int j = k + 1; j <= m; j++)
        {
            s = M[k][j] * x[j]; // формула (4)
            d = d + s; // формула (4)
        }
        x[k] = (B[k] - d) / M[k][k]; // формула (4)
    }
    cout << "Корни системы: " << endl;

    double k_A = x[1], k_b = x[2], k_c = x[3], k_a = exp(k_A), Y_reg2[DIM2], SSE2 = 0, MAE2 = 0;
    
    cout << endl << "a: " << k_a << " b: "<< k_b << " c: " << k_c << endl;
    cout << endl << "Подобранная функциональная зависимость заданного вида: y = (" << k_a << "*x^" << k_b << ")*(e^" << k_c <<"*x)"<< endl;


    for (int i = 0; i < DIM2; i++)
        Y_reg2[i] = (k_a * pow(arr[0][i], k_b)) * exp(k_c * arr[0][i]);
    

    for (int i = 0; i < DIM2; i++)
        SSE2 += pow((arr[1][i] - Y_reg2[i]), 2);
    

    cout << endl << "Суммарная квадратичная ошибка для подобранной функциональной зависимости SSE: " << SSE2 << endl;

    MAE2 = SSE2 / DIM2;

    cout << endl << "Средняя ошибка для подобранной функциональной зависимости MAE: " << MAE2 << endl << endl;

    //Формирование массивов для изображения подобранной функциональной зависимости на графике.
    float y3[DIM2];
    for (int i = 0; i < DIM2; i++)
        y3[i] = (k_A* pow(arr[0][i], k_b)) * exp(k_c * arr[0][i]);

    ofstream f;			
    f.open("data.txt");	
    for (int i=0; i < DIM2; i++) {
        f << arr[0][i] << " ";	
        f << y2[i] << " ";
        f << y3[i] << " ";
        f << arr[1][i] << "\n"; 
    }
    f.close();
}

#### Вывод

![image](https://github.com/pococonut/summer_practice/assets/114181600/395e1c00-355f-463e-ad46-87ef1507ea38)

![image](https://github.com/pococonut/summer_practice/assets/114181600/e4795a28-bb69-43e9-902e-74d3ce24bfac)

![image](https://github.com/pococonut/summer_practice/assets/114181600/a63a4947-604c-4d8b-8195-fabd5334fec6)


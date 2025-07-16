#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

enum Plot_var{
    dIL1 = 0, dUC1, dUC2, dUCb1, dUCb2,
    UL1, UR1, URb1, URu1, UId1,  URb2, URu2, UId2,
    IE1, IE2, IC1, IC2, ICb1, ICb2,
    IL1, IR1, IRb1, IRu1, IId1,  IRb2, IRu2, IId2,
    UE1, UE2, UC1, UC2, UCb1, UCb2,
};

#define VARIABLE IL1   // Переменная, которая выводится на графике

constexpr char RESULT_FILE[] = "result.dat";

#define T_end 1e-3   // Время расчета
#define ACR 1e-3     // Точность
#define EPS_MIN 1e-3 // Минимальное приращение за шаг
#define EPS_MAX 7e-2 // Максимальное приращение за шаг
#define MAX_ITERATION_NEWTON 7  // Максимальное число итераций метода Ньютона

#define SST 1e-7  // Начальный шаг по времени
#define SMN 1e-12    // Минимальный шаг по времени
#define SMX 1e-4     // Максимальный шаг по времени

int N_vd = 0;         // Число ветвей дерева
int N_h = 0;          // Число хорд
int N_ps = 0;         // Число переменных состояния

int id_Uh, id_Ih, id_Uvd, id_Ivd = 0; 

void init_id(){
    id_Uh = N_ps;                  // Номер первого напряжения хорд
    id_Ih = N_ps + N_h + N_vd;     // Номер первого тока хорд
    id_Ivd = N_ps + N_h;            // Номер первого тока ветвей дерева
    id_Uvd = N_ps + 2*N_h + N_vd;   // Номер первого напряжения ветвей дерева  
}                                    // ...........для вычисления индексов всех базисных переменных

Vector create_vector(int size) {
    Vector res;
    res.resize(size, 0);
    return res;
}

Matrix create_matrix(int rows, int cols) {
    Matrix res;
    res.resize(rows);
    for (int i = 0; i < rows; ++i) {
        res[i] = create_vector(cols);
    }
    return res;
}

void zero(Vector& v) {
    fill(v.begin(), v.end(), 0);
}

void zero(Matrix& m) {
    for (auto & i : m) {
        zero(i);
    }
}

// Метод Гаусса для решения СЛАУ
int gauss(Matrix& m, Vector& v) {
    for (int k = 0; k < m.size(); ++k) {
        if (fabs(m[k][k]) < 1e-17) {
            return 1;
        }
        double diagonal = m[k][k];
        // divide this row by diagonal element
        for (int i = k; i < m.size(); ++i) {
            m[k][i] /= diagonal;
        }
        v[k] /= diagonal;
        for (int i = k + 1; i < m.size(); ++i) {
            double elem = m[i][k];
            for (int j = k; j < m.size(); ++j) {
                m[i][j] -= elem * m[k][j];
            }
            v[i] -= elem * v[k];
        }
    }
    for (int i = m.size() - 2; i >= 0; --i) {
        for (int j = i + 1; j < m.size(); ++j) {
            v[i] -= m[i][j] * v[j];
        }
    }
    return 0;
}

// Класс источника напряжения
class E {
public:
    int i, j;
    double val;  // Величина напряжения
    int id;      // Индекс базисной переменной
    int d_type;  // 0 если хорда, 1 если ветвь дерева
    E(int _i, int _j, double _val, int _d_type) {
        i = _i;
        j = _j;
        val = _val;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
    }
    double get_dif_component_UE() { // Производная от невязки для Е
        return 1;
    }

    double get_component(Vector& X) {   // Невязка для обычного источника напряжения
        return X[id+N_ps + d_type*(2*N_h+N_vd)] - val;
    }
};

// Класс синусоидального источника напряжения
class E_sin {
public:
    int i, j;
    double A;  // Величина напряжения
    double P;
    int id;      // Индекс базисной переменной
    int d_type;  // 0 если хорда, 1 если ветвь дерева
    E_sin(int _i, int _j, double _A, double _P, int _d_type) {
        i = _i;
        j = _j;
        A = _A;
        P = _P;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
    }
    double get_dif_component_UE() { // Производная от невязки для Е
        return 1;
    }

    double get_component(Vector& X, double t) {  // Невязка для синусоидального источника
        return X[id+N_ps + d_type*(2*N_h+N_vd)] - A * sin(2 * M_PI / P * t);
    }
};

// Класс источника тока
class I {
public:
    int i, j;
    double val;
    int id;
    int d_type;
    I(int _i, int _j, double _val, int _d_type) {
        i = _i;
        j = _j;
        val = _val;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
    }
    double get_dif_component_II() {
        return 1;
    }

    double get_component(Vector& X) {
        return X[id+id_Ih - d_type*N_vd] - val;
    }
};

// Класс источника тока для диода
class Id {
public:
    int i, j;
    int id;
    int d_type;
    double It, MFt;
    Id(int _i, int _j, double _It, double _MFt, int _d_type) {
        i = _i;
        j = _j;
        It = _It;
        MFt = _MFt;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
    }
    double get_dif_component_IId() {
        return 1;
    }
    double get_dif_component_UId(Vector& X) {
        if (d_type)
            return -It / MFt * exp( X[id+2*N_h+N_vd+N_ps]/ MFt);
        return -It / MFt * exp( X[id+N_ps]/ MFt);
    }
    double get_component(Vector& X) {
        return X[id+id_Ih - d_type*N_vd] - It * (exp(X[id+N_ps + d_type*(2*N_h+N_vd)] / MFt) - 1);
    }
};

// Класс конденсатора
class C {
public:
    int i, j;
    double val;
    int id;
    int id_ps;  // Индекс производной переменной состояния для конденсатора
    int d_type;
    C(int _i, int _j, double _val, int _d_type) {
        i = _i;
        j = _j;
        val = _val;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
        id_ps = N_ps;
        N_ps++;
    } 
    double get_dif_component_dUC() {  // Производная от невязки по dUC/dt
        return -val;
    }
    double get_dif_component_IC() {  // Производная от невязки по IC
        return 1;
    }
    double get_component(Vector& X) {  // Невязка для конденсатора (компонентное уравнение)
        return X[id+id_Ih - d_type*N_vd] - val* X[id_ps];
    }
};

// Класс катушки
class L {
public:
    int i, j;
    double val;
    int id;
    int id_ps;
    int d_type;
    L(int _i, int _j, double _val, int _d_type) {
        i = _i;
        j = _j;
        val = _val;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
        id_ps = N_ps;
        N_ps++;
    }
    double get_dif_component_dIL() {
        return -val;
    }
    double get_dif_component_UL() {
        return 1;
    }
    double get_component(Vector& X) {
        return X[id+N_ps + d_type*(2*N_h+N_vd)] - val* X[id_ps];
    }
};

// Класс резистора
class R {
public:
    int i, j;
    double val;
    int id;
    int d_type;
    R(int _i, int _j, double _val, int _d_type) {
        i = _i;
        j = _j;
        val = _val;
        d_type = _d_type;

        if (_d_type) {
            id = N_vd;
            N_vd++;}
        else{
            id = N_h;
            N_h++;
        }
    }
    double get_dif_component_IR() {
        if (d_type)
            return -1;
        return val;
    }
    double get_dif_component_UR() {
        if (d_type)
            return 1/val;
        return -1;
    }
    double get_component(Vector& X) {
        if (d_type)
            return X[id+id_Uvd]/val - X[id+id_Ivd];
        return X[id+id_Ih]* val - X[id + id_Uh];
    }
};

double find_max_elem(const Vector& v) {
    double m = -1;
    for (auto& e : v) {
        if (fabs(e) > m) {
            m = fabs(e);
        }
    }
    return m;
}

// Класс схемы
class Scheme {
public:
    vector<C> Cs;  // Вектора для каждого типа элементов
    vector<L> Ls;
    vector<R> Rs;
    vector<Id> Ids;
    vector<E> Es;
    vector<E_sin> Esins;
    vector<I> Is;

    Matrix M_matrix(int N_nodes){  // Функция построения М-матрицы контуров и сечений
        Matrix graph_m = create_matrix(N_h+N_vd, N_nodes+1);  // Матрица соответствия узлов и элементов
        Matrix M = create_matrix(N_h, N_vd);  // М-матрица
        zero(graph_m);

        // Заполнение матрицы соответствия узлов и элементов
        for (auto& e : Es) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
         for (auto& e : Esins) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
        for (auto& e : Is) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
        for (auto& e : Ids) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
        for (auto& e : Rs) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
        for (auto& e : Ls) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }
        for (auto& e : Cs) {
            graph_m[e.id+e.d_type*N_h][e.i] = 1; graph_m[e.id+e.d_type*N_h][e.j]=-1;
        }

        Matrix graph_m_1 = graph_m;
        int i_prev = 0;

        // Поиск циклов в графе и заполнение М-матрицы
        for (int i=0; i<N_h; i++){
            int j_prev = i;
            if (i != i_prev){ graph_m_1 = graph_m;}
            i_prev = i;

            int i_end = find(begin(graph_m_1[i]), end(graph_m_1[i]), 1) - begin(graph_m_1[i]);
            int i_start = find(begin(graph_m_1[i]), end(graph_m_1[i]), -1) - begin(graph_m_1[i]);

            int iter = 0;
            
            zero(M[i]);
            // Для каждой хорды ищем ветви дерева с которыми образуется цикл. 
            // В М-матрицу заносится -1 или 1 для ветви дерева, в зависимости от направления
            while(i_start != i_end){
                iter++;
                for (int j=N_h; j<graph_m_1.size(); j++){
                    if (graph_m_1[j][i_start] != 0) {
                        if (j != j_prev){
                            j_prev = j;
                            M[i][j-N_h] = graph_m_1[j][i_start];
                            i_start = find(begin(graph_m_1[j]), end(graph_m_1[j]), -graph_m_1[j][i_start]) - begin(graph_m_1[j]);
                            break;
                        }
                    }
                    if ((j == graph_m_1.size()-1) || (iter > N_nodes)) { // Если цикл не найден, переходим к поиску среди других ветвей дерева
                        if (j != j_prev){
                            zero(graph_m_1[j_prev]);
                            i = i_prev-1;
                            i_start = i_end;
                            break;
                        }
                    }
                }
            }
        }

        return M;
    }

    void fill_matrix(Matrix& Y, Vector& X_cur, Matrix& M, Matrix& M_T, double dt){
        // Заполнение матрицы производными от разностных схем
        for (auto& e : Cs) {
            int ind = 0;
            if (e.d_type) ind = id_Uvd;
            else ind = id_Uh;
            Y[e.id_ps][e.id_ps] = 1;
            Y[e.id_ps][e.id+ind] = -1/dt;
        }
        for (auto& e : Ls) {
            int ind = 0;
            if (e.d_type) ind = id_Ivd;
            else ind = id_Ih;
            Y[e.id_ps][e.id_ps] = 1;
            Y[e.id_ps][e.id+ind] = -1/dt;
        }

        // Заполнение матрицы производными от топологических уравнений
        for (int i=N_ps; i<N_h+N_vd+N_ps; i++){
            Y[i][i] = 1.0;
        }
        for (int i=id_Uh; i<id_Ivd; i++){
            for (int j=id_Uvd; j<Y.size(); j++){
                Y[i][j] = M[i-id_Uh][j-id_Uvd];
            }
        }
        for (int i=id_Ivd; i<id_Ih; i++){
            for (int j=id_Ih; j<id_Uvd; j++){
                Y[i][j] = -1*M_T[i-id_Ivd][j-id_Ih];
            }
        }

        // Заполнение матрицы производными от компонентных уравнений
        for (auto& e : Cs) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Ih - e.d_type*N_vd] = e.get_dif_component_IC();
            Y[e.id+id_Ih + e.d_type*N_h][e.id_ps] = e.get_dif_component_dUC();
        }
        for (auto& e : Ls) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Uh + e.d_type*(2*N_h+N_vd)] = e.get_dif_component_UL();
            Y[e.id+id_Ih + e.d_type*N_h][e.id_ps] = e.get_dif_component_dIL();
        }
        for (auto& e : Rs) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Uh + e.d_type*(2*N_h+N_vd)] = e.get_dif_component_UR();
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Ih - e.d_type*N_vd] = e.get_dif_component_IR();
        }
        for (auto& e : Ids) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Uh + e.d_type*(2*N_h+N_vd)] = e.get_dif_component_UId(X_cur);
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Ih - e.d_type*N_vd] = e.get_dif_component_IId();
        }
        for (auto& e : Is) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Ih - e.d_type*N_vd] = e.get_dif_component_II();
        }
        for (auto& e : Es) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Uh + e.d_type*(2*N_h+N_vd)] = e.get_dif_component_UE();
        }
        for (auto& e : Esins) {
            Y[e.id+id_Ih + e.d_type*N_h][e.id+id_Uh + e.d_type*(2*N_h+N_vd)] = e.get_dif_component_UE();
        }   
    }

    void fill_vector(Vector& b, Vector& X_cur, Vector& X_prev, Matrix& M, Matrix& M_T, double dt, double t){

        // Заполнение вектора невязок разностными схемами
        for (auto& e : Cs) {
            int UC;
            if (e.d_type) UC = e.id + id_Uvd;
            else UC = e.id + id_Uh;
            b[e.id_ps] = X_cur[e.id_ps] - (X_cur[UC]-X_prev[UC])/dt;
        }
        for (auto& e : Ls) {
            int IL;
            if (e.d_type) IL = e.id + id_Ivd;
            else IL = e.id + id_Ih;
            b[e.id_ps] = X_cur[e.id_ps] - (X_cur[IL]-X_prev[IL])/dt;
        }

        // Заполнение вектора невязок топологическими уравнениями
        for (int i=0; i < M.size(); i++){
            b[i+id_Uh] = X_cur[i+id_Uh];
            for (int j=0; j < M[0].size(); j++){
                b[i+id_Uh] += M[i][j]*X_cur[j+id_Uvd];
            }
        }
        for (int i=0; i < M_T.size(); i++){
            b[i+id_Ivd] = X_cur[i+id_Ivd];
            for (int j=0; j < M_T[0].size(); j++){
                b[i+id_Ivd] -= M_T[i][j]*X_cur[j+id_Ih];
            }
        }

        // Заполнение вектора невязок компонентными уравнениями
        for (auto& e : Cs) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Ls) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Rs) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Ids) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Is) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Es) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur);
        }
        for (auto& e : Esins) {
            b[e.id+id_Ih + e.d_type*N_h] = e.get_component(X_cur, t);
        }
    }

    void calculate(double T, double epsilon, double eps_min, double eps_max, int N_nodes) {
        // Файлы для записи результатов
        ofstream result_file(RESULT_FILE);
        result_file.setf(ios_base::right);
        init_id();

        Matrix M = M_matrix(N_nodes);  // Заполнение М-матрицы
        Matrix M_T(M[0].size());    // Транспонирование М-матрицы
        for (int i = 0; i < M[0].size(); i++)
        {
            for (int j = 0; j < M.size(); j++)
            {
                M_T[i].push_back(M[j][i]);
            }
        }

        cout << "M-matrix" << endl;  // Вывод М-матрицы
        for (int i=0; i<M.size(); i++){
            for (int j=0; j<M[0].size(); j++){
                cout << M[i][j] <<" ";}
            cout << endl;
        }
        cout << endl;

        int count = N_ps + 2*(N_vd+N_h);  // Общее число уравнений
        Matrix Y = create_matrix(count, count);  // Матрица Якоби
        Vector b = create_vector(count);  // Вектор невязок

        Vector X_cur = create_vector(count);  // Вектор базисных переменных
        Vector X_prev = create_vector(count);  // Значения базисных переменных на предыдущем шаге
        Vector X_prev_prev = create_vector(count);  // Значения базисных переменных 2 шага назад

        double step_count = 0;  // Сумма шагов
        double iter_count = 0;  // Сумма итераций
        int step = 0;  // Счетчик шагов

        double dt_prev = 0;  
        double dt = SST;  // Шаг по времени

        // Цикл итераций по времени
        double t = dt;
        while (t <= T) {
            if (step % 100 == 0 ) cout << "t=" << t << endl;

            bool solved = false; 
            int n = 0;
            // Цикл итераций метода Ньютона
            while (!solved) {
                zero(Y);
                fill_matrix(Y, X_cur, M, M_T, dt);  // Заполнение матрицы и вектора
                fill_vector(b, X_cur, X_prev, M, M_T, dt, t);

                for (auto& e : b) {
                    e *= -1;
                }

                //Решение СЛАУ
                gauss(Y, b);

                // Обновление значений базисных переменных
                for (int i = 0; i < X_cur.size(); ++i) {
                    X_cur[i] += b[i];
                }

                double dX_max = find_max_elem(b);
                if (dX_max < epsilon) solved = true;  // Индикатор сходимости метода Ньютона
                n++;

                if (!solved){  // Если метод не сошелся (точность не достигнута)
                    if (n>7){   // Если число итераций больше максимума
                        n=0; 
                        dt *= 0.5;   // Уменьшение шага
                        X_cur = X_prev;  // Отбрасываем результаты текущего шага
                    
                        if (dt < SMN) {
                            throw domain_error("Решение не сходится, dt < dt_min");
                        }
                    }
                }
            }
            double cur_delta = 0.0;
            // Расчет текущего приращения
            for (int i = id_Uh; i < X_cur.size(); i++) {
                double tmp = abs(0.5 * dt * ( (X_cur.at(i) - X_prev.at(i) ) / dt - (X_prev.at(i) - X_prev_prev.at(i)) / dt_prev));
                cur_delta = (tmp > cur_delta) ? tmp : cur_delta;
            }
            // Если интегрирование неудовлетворительно по точности
            if (cur_delta > eps_max && dt>SMN) {
                if (dt < SMN) {
                    throw domain_error("Решение не сходится, cur_delta > eps_max");
                }
                // Уменьшаем шаг по времени и отбрасываем результаты текущего шага
                dt *= 0.5;
                X_cur = X_prev;
            } 
            else {
                // Сохранение значений с предыдущего шага
                X_prev_prev = X_prev;
                X_prev = X_cur;
                dt_prev = dt;
                // Вывод значения безисных переменных на текущем временном шаге
                result_file << setw(15) << t;
                for (auto& e : X_cur)
                    result_file << setw(15) << e;
                result_file << endl;

                t += dt;
                step++;
                step_count += dt;
                iter_count+=n;

                if (cur_delta < eps_min || dt<SMN) {  // Если приращение слишком мало
                    dt *= 2;
                } 
            }
        }
        cout << "\nN_steps: " << step << endl;
        cout << "Average step: " << step_count/step << endl;
        cout << "Average number of iterations: " << iter_count/step << endl;
    }
};

// Вывод графика
void plot() {
    ofstream gnuscript_file("plt.gnuplot");
    gnuscript_file << "set key right bottom\n";
    gnuscript_file << "set terminal wxt size 900,600\n";
    gnuscript_file << "set grid\n";
    gnuscript_file << "set xrange[" << 0 << ':' << T_end << ']' << endl;
    gnuscript_file << "plot '" << RESULT_FILE << "' using 1:" << VARIABLE+2 <<" with l lc rgb '#000000' lw 1.5 title ' '" << endl;
    gnuscript_file << "pause -1" << endl;
    system((string("gnuplot ") + string("plt.gnuplot")).c_str());
}

int main() {
    Scheme scheme;

    int N_nodes = 6; // Число узлов в схеме (не включая базу) 

    // Включаем элементы схемы: вызываем функцию для определенного типа элемента:
    // 1й аргумент -- узел, из которого выходит ток, 2й аргумент -- узел, в который ток втекает,  (0-базовый узел)
    // 3й -- величина напряжения/емкости/индуктивности/сопротивления, последний аргумент -- 1 если элемент входит в ветви дерева графа, иначе 0
    scheme.Ls.push_back(L(5, 6, 2.53e-4, 0));        // L1
    scheme.Rs.push_back(R(5, 0, 1000, 0));           // R1
    scheme.Rs.push_back(R(2, 4, 20, 0));             // Rb1
    scheme.Rs.push_back(R(4, 5, 1000000, 0));        // Ru1
    scheme.Ids.push_back(Id(4, 5, 1e-12, 0.026, 0)); // Id1
    scheme.Rs.push_back(R(0, 3, 20, 0));             // Rb2
    scheme.Rs.push_back(R(3, 2, 1000000, 0));        // Ru2
    scheme.Ids.push_back(Id(3, 2, 1e-12, 0.026, 0)); // Id2
    
    scheme.Esins.push_back(E_sin(1, 0, 10, 1e-4, 1));      // E1
    scheme.Es.push_back(E(6, 0, 5, 1));            // E2
    scheme.Cs.push_back(C(1, 2, 1e-6, 1));        // C1
    scheme.Cs.push_back(C(5, 0, 1e-6, 1));        // C2
    scheme.Cs.push_back(C(4, 5, 2e-12, 1));       // Cb1
    scheme.Cs.push_back(C(3, 2, 2e-12, 1));       // Cb2

    scheme.calculate(T_end, ACR, EPS_MIN, EPS_MAX, N_nodes);

    plot();

    return 0;
}

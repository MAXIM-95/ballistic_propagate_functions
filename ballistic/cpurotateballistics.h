#pragma once

#include "math/ballist/cpuballistics.h"

namespace RMath {

using namespace RMatrix;
using namespace RModels::RAtmosphere::RAtmInterfaces;

namespace RBallistics {

/// \brief Дополнительные параметры элемента для расчёта сил и моментов
struct CPUElemParams
{
protected:
  T    wx,                ///< параметры вращения блока
       wy,                ///< параметры вращения блока
       wz;                ///< параметры вращения блока
  //... расчитываемые параметры ...
  T    _gravityCenter,    ///< центр масс
       _pressureCenter;   ///< центр давления

  /*!
   * \brief Расстояние от носика до центра тяжести вдоль продольной оси
   * \return Расстояние от носика до центра тяжести вдоль продольной оси
   */
  T centerOfGravityShift();

public:

  T    weight,            ///< масса
       length,            ///< длина
       base_diameter,     ///< диаметр основания
       midel,             ///< площадь Миделя
       nose_radius,       ///< радиус "носика"
       aerodynamic_koef,  ///< аэродин. коэф. вдоль продольной оси скоростной СК
                          ///< при нулевом угле атаки
       polara_koef,       ///< коэф. поляры
       lifting_force,     ///< производная коэф. подъёмной силы по углу атаки в скоростной СК
       inert_moment_x,    ///< главный момент инерции
       inert_moment_y,    ///< главный момент инерции
       inert_moment_z,    ///< главный момент инерции
       pitch_momen_koef,  ///< коэф. момента тангажа

       attack_angle;      ///< угол атаки на входе в атмосферу

  /*!
  * \param iwx Скорость закрутки элемента
  * \param iwy Скорость закрутки элемента
  * \param iwz Скорость закрутки элемента
  */
    CPUElemParams(const T &iwx, const T &iwy, const T &iwz);
    inline CPUElemParams(const CPUElemParams &other){ *this = other; }
    CPUElemParams & operator=(const CPUElemParams &other);

  /*!
  * \brief Расчёт начальных углов на входе в атмосферу
  * \param vx вектор скорости
  * \param vy вектор скорости
  * \param vz вектор скорости
  * \param pitch Тангаж
  * \param raw Рыскание
  * \param banking Крен
  * \param owx Вектор угловой скорости вращения объекта
  * \param owy Вектор угловой скорости вращения объекта
  * \param owz Вектор угловой скорости вращения объекта
  * \return Начальные углы на входе в атмосферу
  */
  bool calculateInitialAngles(const T &vx, const T &vy, const T &vz,
                              T &pitch, T &raw, T &banking, T &owx, T &owy, T &owz);

};


class CPUMPSKRotatePropagator: public CPUBasePropagator{
private:
  CPUElemParams      *_params;
  T                   _centerx,   ///< координаты центра МПСК в ПГСК
                      _centery,   ///< координаты центра МПСК в ПГСК
                      _centerz;   ///< координаты центра МПСК в ПГСК
  T                   _sinfi,     ///< функции широты точки стояния
                      _cosfi;     ///< функции широты точки стояния
  T                   _maxheight; ///< пороговое значение высоты атмосферы
  T                   _currency;  ///< точность расчёта пороговой высоты

  /*!
   * \brief Создание и удаление калькулятора траекторий в МПСК (Моисеевский)
   * \param ycsk Координата y объекта в ПГСК
   * \param zcsk Координата z объекта в ПГСК
   * \param dist расстояние от объекта до центра земли
   * \return Высота объекта над поверхностью эллипсоида
   */
  CPUMPSKRotatePropagator(const CPUMPSKRotatePropagator&);
  CPUMPSKRotatePropagator &operator=(const CPUMPSKRotatePropagator&);

protected:
  /*!
   * \brief Расчет высоты объекта над поверхностью эллипсоида по координатам в МПСК
   * \param ycsk Координата y объекта в ПГСК
   * \param zcsk Координата z объекта в ПГСК
   * \param dist расстояние от объекта до центра земли
   * \return Высота объекта над поверхностью эллипсоида
   */
  T _heightInMPSK(const T &ycsk, const T &zcsk, const T &dist);
  /*!
   * \brief Расчет полного ускорения баллистического объекта
   * \param x вектор положения
   * \param y вектор положения
   * \param z вектор положения
   * \param vx вектор скорости
   * \param vy вектор скорости
   * \param vz вектор скорости
   * \param pitch Тангаж
   * \param raw Рыскание
   * \param banking Крен
   * \param dvx вектор ускорений
   * \param dvy вектор ускорений
   * \param dvz вектор ускорений
   * \param gamma Баллистический коэффициент
   * \param gamma1 Баллистический коэффициент
   * \param inatm Признак пребывания в атмосфере
   * \return Флаг выполнения расчётов
   */
  bool accelerations(const T &x, const T &y, const T &z,
                     const T &vx, const T &vy, const T &vz,
                     T &dvx, T &dvy, T &dvz,
                     const T &pitch, const T &raw, const T &banking,
                     T &vpitch, T &vraw, T &vbanking,
                     T &wx, T &wy, T &wz,
                     T &dwx, T &dwy, T &dwz, bool inatm, T &hght,
                     const T &atmorange=NAN);

  /*!
   * \brief Преобразование МПСК ==> Связанная СК
   * \param pitch Тангаж
   * \param raw Рыскание
   * \param banking Крен
   * \return
   */
  CPUMatrix *MPSK2LSK(const T &pitch, const T &raw, const T &banking);

    /*!
   * \brief Прогноз движения объекта в поле Земли на один шаг (в МПСК)
   * \param t время привязкии исходных координат
   * \param dt Шаг интегрирования
   * \param x вектор положения
   * \param y вектор положения
   * \param z вектор положения
   * \param vx вектор скоростей
   * \param vy вектор скоростей
   * \param vz вектор скоростей
   * \param dvx вектор ускорений
   * \param dvy вектор ускорений
   * \param dvz вектор ускорений
   * \param pitch Тангаж
   * \param raw Рыскание
   * \param banking Крен
   * \param vpitch Скорость изменения тангажа
   * \param vraw Скорость изменения рыскания
   * \param vbanking Скорость изменения крена
   * \param wx Вектор угловой скорости вращения объекта
   * \param wy Вектор угловой скорости вращения объекта
   * \param wz Вектор угловой скорости вращения объекта
   * \param dwx Производная вектора угловой скорости вращения объекта
   * \param dwy Производная вектора угловой скорости вращения объекта
   * \param dwz Производная вектора угловой скорости вращения объекта
   * \param inatm Признак пребывания в атмосфере
   * \param hght Рассчитанное значение высоты
   * \param atmorange Порог учёта плотности атмосферы
   * \return Результат выполнения операции
   */
  PropagateResult __propagate(T &t, const T &dt, T &x, T &y, T &z, T &vx, T &vy, T &vz,
                              T &dvx, T &dvy, T &dvz, T &pitch, T &raw, T &banking,
                              T &vpitch, T &vraw, T &vbanking, T &wx, T &wy, T &wz,
                              T &dwx, T &dwy, T &dwz, bool inatm, T &hght,
                              const T &atmorange=NAN);
public:
  /*!
   * \brief Высота над поверхностью эллипсоида
   * \param x вектор положения
   * \param y вектор положения
   * \param z вектор положения
   * \return
   */
  virtual T heightInBaseSystem(const T &x, const T &y, const T &z);

    /*!
    * \brief Функция расчёта скорости изменения высоты
    * \param x составляющая положения в МПСК
    * \param y составляющая положения в МПСК
    * \param z составляющая положения в МПСК
    * \param vx составляющая скорости в МПСК
    * \param vy составляющая скорости в МПСК
    * \param vz составляющая скорости в МПСК
    */
  T altitude_rate(const T &x, const T &y, const T &z,
                  const T &vx,const T &vy,const T &vz);

  /*!
   * \brief CPUMPSKPropagator
   * \param ib Широта стояния МПСК, [рад]
   * \param ih Высота стояния МПСК, [м]
   * \param imaxheight Пороговое значение высоты атмосферы, [м]
   * \param iprms Дополнительные параметры элемента для расчёта сил и моментов
   * \param atm_type Тип атмосферы (SA_UNUSED - не используется)
   * \param iseason Сезон для расчёта плотности атмосферы
   * \param icurr Точность достижения высоты
   * \param dfront   Величина рубежа
   * \param frx Координата X центра рубежа
   * \param fry Координата Y центра рубежа
   * \param frz Координата Z центра рубежа
   */
  explicit CPUMPSKRotatePropagator(const T& ib, const T& ih,
                                   const T& imaxheight,
                                   IStaticAtmosphere::ISATypes atm_type,
                                   const TEarthDatum &iearth,
                                   rint32 iseason=IStaticAtmosphere::S_SUMMER,
                                   T icurr = 1.e-6, const T &dfront=NAN,
                                   const T &frx=NAN, const T &fry=NAN,
                                   const T &frz=NAN);
  inline virtual ~CPUMPSKRotatePropagator(){}

    /*!
   * \brief Прогноз движения объекта в поле Земли на один шаг (в МПСК)
   * \param iprms Параметры объекта, описывающие вращение вокруг центра масс
   * \param t время привязкии исходных координат
   * \param dt Шаг интегрирования
   * \param x вектор положения
   * \param y вектор положения
   * \param z вектор положения
   * \param vx вектор скоростей
   * \param vy вектор скоростей
   * \param vz вектор скоростей
   * \param dvx вектор ускорений
   * \param dvy вектор ускорений
   * \param dvz вектор ускорений
   * \param pitch Тангаж
   * \param raw Рыскание
   * \param banking Крен
   * \param vpitch Скорость изменения тангажа
   * \param vraw Скорость изменения рыскания
   * \param vbanking Скорость изменения крена
   * \param wx Вектор угловой скорости вращения объекта
   * \param wy Вектор угловой скорости вращения объекта
   * \param wz Вектор угловой скорости вращения объекта
   * \param dwx Производная вектора угловой скорости вращения объекта
   * \param dwy Производная вектора угловой скорости вращения объекта
   * \param dwz Производная вектора угловой скорости вращения объекта
   * \param inatm Признак пребывания в атмосфере
   * \param hght Рассчитанное значение высоты
   * \param atmorange Порог учёта плотности атмосферы
   * \return Результат выполнения операции
   */
  virtual PropagateResult propagate(CPUElemParams *iprms, T &t, const T &dt,
                                    T &x, T &y, T &z,
                                    T &vx, T &vy, T &vz,
                                    T &dvx, T &dvy, T &dvz,
                                    T &pitch, T &raw, T &banking,
                                    T &vpitch, T &vraw, T &vbanking,
                                    T &wx, T &wy, T &wz,
                                    T &dwx, T &dwy, T &dwz, bool &inatm,
                                    T &hght, const T &atmorange=NAN);

};




} //namespace RBallistics

} //namespace RMath


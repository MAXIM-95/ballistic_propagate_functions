#pragma once

#include "math/cpumath.h"
#include "math/geo/cpugeomath.h"
#include "math/matrix/cpumatrix.h"
#include "math/ballist/cpuballistics.h"
#include "models/datums/cpugeodatums.h"
#include "models/atmo/staticatmosphere.h"

namespace RMath {

using namespace RModels::RDatums;
using namespace RModels::RAtmosphere::RAtmInterfaces;

namespace RBallistics {


class CPU45Propagator: public CPUBasePropagator{
private:
  CPU45Propagator(const CPU45Propagator&);
  CPU45Propagator &operator=(const CPU45Propagator&);

protected:
  /*!
  * \brief Коэффициент №1
  * \param x вектор положения
  * \param y вектор положения
  * \param z вектор положения
  * \param vx вектор скорости
  * \param vy вектор скорости
  * \param vz вектор скорости
  * \param dt Шаг интегрирования
  * \param atmorange Порог учёта плотности атмосферы
  * \return Плотность атмосферы на высоте height
  */
  T K1(const T &x, const T &y, const T &z,
       const T &Vx, const T &Vy, const T &Vz, const T &dt, const T &bk,
       const T &atmorange=NAN);
  /*!
  * \brief Коэффициент №2
  * \param x вектор положения
  * \param y вектор положения
  * \param z вектор положения
  * \param vx вектор скорости
  * \param vy вектор скорости
  * \param vz вектор скорости
  * \param dt Шаг интегрирования
  * \param atmorange Порог учёта плотности атмосферы
  * \return Плотность атмосферы на высоте height
  */
  T K2(const T &x, const T &y, const T &z,
       const T &Vx, const T &Vy, const T &Vz, const T &dt, const T &bk,
       const T &atmorange=NAN);
  /*!
  * \brief Коэффициент №3
  * \param x вектор положения
  * \param y вектор положения
  * \param z вектор положения
  * \param vx вектор скорости
  * \param vy вектор скорости
  * \param vz вектор скорости
  * \param dt Шаг интегрирования
  * \param atmorange Порог учёта плотности атмосферы
  * \return Плотность атмосферы на высоте height
  */
  T K3(const T &x, const T &y, const T &z,
       const T &Vx, const T &Vy, const T &Vz, const T &dt, const T &bk,
       const T &atmorange=NAN);

public:
  /*!
   * \brief CPUMPSKPropagator
   * \details Cоздание и удаление калькулятора траекторий в ПГСК
   * \param ib Широта стояния МПСК, [рад]
   * \param ih Высота стояния МПСК, [м]
   * \param ihstab Высота стабилизации элемента, [м]
   * \param atm_type Тип атмосферы (SA_UNUSED - не используется)
   * \param iseason Сезон для расчёта плотности атмосферы
   * \param iearth Параметры Земли
   * \param dfront   Величина рубежа
   * \param frx Координата X центра рубежа
   * \param fry Координата Y центра рубежа
   * \param frz Координата Z центра рубежа
   */
  explicit CPU45Propagator(IStaticAtmosphere::ISATypes atm_type, rint32 iseason,
                           const TEarthDatum &iearth, const T &dfront=NAN,
                           const T &frx=NAN, const T &fry=NAN, const T &frz=NAN);
  inline virtual ~CPU45Propagator(){}

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
   * \brief Прогноз движения объекта в поле Земли на один шаг (в МПСК)
   * \param t время привязки исходных координат
   * \param dt Шаг интегрирования
   * \param x вектор положения
   * \param y вектор положения
   * \param z вектор положения
   * \param vx вектор ускорений
   * \param vy вектор ускорений
   * \param vz вектор ускорений
   * \param dvx вектор ускорений
   * \param dvy вектор ускорений
   * \param dvz вектор ускорений
   * \param gamma Баллистический коэффициент
   * \param gamma1 Баллистический коэффициент
   * \param atmorange Порог учёта плотности атмосферы
   * \return Результат выполнения операции
   */
  virtual PropagateResult propagate(T &t, const T &dt,
                         T &x, T &y, T &z,
                         T &vx, T &vy, T &vz,
                         T &dvx, T &dvy, T &dvz,
                         const T &gamma, const T &gamma1,
                         const T &atmorange=NAN);


  /*!
   * \brief Прогноз движения объекта в поле Земли до достижения высоты (в МПСК)
   * \param t время привязкии исходных координат в сек
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
   * \param height пороговое значение высоты, до кот. проводится моделирование
   * \param idir участок траектории, для кот. задано пороговое значение
   * \param maxflighttime макс. длительность прогноза в сек
   * \param atmorange Порог учёта плотности атмосферы
   * \return Результат выполнения операции
   */
    virtual PropagateResult propagateOnHeight(T &t, const T &dt,
                                              T &x, T &y, T &z,
                                              T &vx, T &vy, T &vz,
                                              T &dvx, T &dvy, T &dvz,
                                              const T &gamma, const T &gamma1,
                                              T &height, FlightDirection idir,
                                              const T &maxflighttime,
                                              const T &atmorange=NAN);


};



} //namespace RBallistics

} //namespace RMath


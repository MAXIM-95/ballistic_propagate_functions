#pragma once

#include "math/cpumath.h"
#include "math/geo/cpugeomath.h"
#include "math/matrix/cpumatrix.h"
#include "models/datums/cpugeodatums.h"
#include "math/ballist/cpuballistics.h"
#include "models/atmo/staticatmosphere.h"

namespace RMath {

using namespace RModels::RDatums;
using namespace RModels::RAtmosphere::RAtmInterfaces;

namespace RBallistics {

class CPUPGSKPropagator: public CPUBasePropagator{
private:
  T                   _rostab;  ///< Плотность атмосферы на высоте стабилизации

  CPUPGSKPropagator(const CPUPGSKPropagator&);
  CPUPGSKPropagator &operator=(const CPUPGSKPropagator&);

protected:
  /*!
   * \brief Расчет полного ускорения баллистического объекта
   * \param x Вектор положения
   * \param y Вектор положения
   * \param z Вектор положения
   * \param vx Вектор скорости
   * \param vy Вектор скорости
   * \param vz Вектор скорости
   * \param dvx Вектор ускорений
   * \param dvy Вектор ускорений
   * \param dvz Вектор ускорений
   * \param gamma Баллистический коэффициент
   * \param gamma1 Баллистический коэффициент
   * \param atmorange Порог учёта плотности атмосферы
   * \return Флаг выполнения расчётов
   */
  bool accelerations(const T &x, const T &y, const T &z,
                     const T &vx, const T &vy, const T &vz, T &dvx, T &dvy, T &dvz,
                     const T &gamma, const T &gamma1, const T &atmorange=NAN);

public:
  /*!
   * \brief Создание и удаление калькулятора траекторий в МПСК
   * \param ib Широта стояния МПСК, [рад]
   * \param ih Высота стояния МПСК, [м]
   * \param ihstab Высота стабилизации элемента, [м]
   * \param atm_type Тип атмосферы (SA_UNUSED - не используется)
   * \param iseason Сезон для расчёта плотности атмосферы
   * \param dfront   Величина рубежа
   * \param frx Координата X центра рубежа
   * \param fry Координата Y центра рубежа
   * \param frz Координата Z центра рубежа
   */
  explicit CPUPGSKPropagator(IStaticAtmosphere::ISATypes atm_type,
                             const TEarthDatum &iearth,
                             const T& ihstab = 80000.,
                             rint32 iseason=IStaticAtmosphere::S_SUMMER,
                             const T &dfront=NAN, const T &frx=NAN,
                             const T &fry=NAN, const T &frz=NAN);
  inline virtual ~CPUPGSKPropagator(){}

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
   * \param t Время привязкии исходных координат
   * \param dt Шаг интегрирования
   * \param x Вектор положения
   * \param y Вектор положения
   * \param z Вектор положения
   * \param vx Вектор скорости
   * \param vy Вектор скорости
   * \param vz Вектор скорости
   * \param dvx Вектор ускорений
   * \param dvy Вектор ускорений
   * \param dvz Вектор ускорений
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
                                   const T &maxflighttime, const T &atmorange=NAN);

};


} //namespace RBallistics

} //namespace RMath


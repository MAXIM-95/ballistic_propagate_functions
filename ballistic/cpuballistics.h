#pragma once

#include "math/cpumath.h"
#include "math/geo/cpugeomath.h"
#include "math/matrix/cpumatrix.h"
#include "models/datums/cpugeodatums.h"
#include "models/atmo/staticatmosphere.h"

namespace RMath {

using namespace RMatrix;
using namespace RModels::RAtmosphere::RAtmInterfaces;

namespace RBallistics {

/*!
 * \brief Интерфейс для методов прогноза
 */
class CPUBasePropagator {
private:
  CPUBasePropagator(const CPUBasePropagator&);
  CPUBasePropagator &operator=(const CPUBasePropagator&);

protected:
  rint32              _season; ///< Сезон для расчёта плотности атмосферы
  TEarthDatum         _datum;  ///< Модель параметров Земли
  IStaticAtmosphere  *_atmodel;///< Используемая модель атмосферы
  T                   _deltaFrontier; ///< отклонение от рубежа
  T                   _frontierX,///центр рубежа
                      _frontierY,///центр рубежа
                      _frontierZ;///центр рубежа

  /*!
  * \brief Расчёт высоты объекта над эллипсоидом
  * \param r Радиус-вектор
  * \param z Координата Z
  * \param datum Параметры Земли
  */
  static T heightInPGSK(const T &r, const T &z,
                        const RModels::RDatums::TEarthDatum &datum);

public:
  ///направление траектории
  enum FlightDirection
  {
     FD_UNDETERMINED=0, ///<восходящий участок траектории
     FD_ASCENDING,      ///<восходящий участок траектории
     FD_DESCENDING      ///<нисходящий участок траектории
  };

  ///результат работы баллистического прогноза
  enum PropagateResult
  {
    PR_OK=0,                     ///< Штатное завершение
    PR_FRONTIER,                 ///< Прерван - выход за указанный рубеж
    PR_OK_WITH_LIMITATIONS,      ///< Прогноз прерван по критическому условию
    PR_INVALID_INPUT_PARAMETERS, ///< Некорректные входные параметры
    PR_MATHERROR,                ///< Прерван - ошибка в вычислениях
    PR_UNIMPLEMENTED             ///< Отсутствует реализация
  };

  static const char * PropagateResultDescription(PropagateResult in);

  /*!
  * \brief Расчёт высоты объекта над эллипсоидом
  * \param x вектор положения
  * \param y вектор положения
  * \param z вектор положения
  * \param datum параметры Земли
  */
  static T heightInPGSK(const T &x, const T &y, const T &z,
                        const RModels::RDatums::TEarthDatum &datum);
///FIXME: ЭТОТ МЕТОД НЕПРАВИЛЬНО ИНТЕРПОЛИРУЕТ!!!
  /*!
   * \brief Интерполяция положения объекта
   * \param t1 Время привязки первого вектора положения и скорости объекта
   * \param x1 Составляющая первого вектора положения объекта
   * \param y1 Составляющая первого вектора положения объекта
   * \param z1 Составляющая первого вектора положения объекта
   * \param vx1 Составляющая первого вектора скорости объекта
   * \param vy1 Составляющая первого вектора скорости объекта
   * \param vz1 Составляющая первого вектора скорости объекта
   * \param t2 Время привязки второго вектора положения и скорости объекта
   * \param x2 Составляющая второго вектора положения объекта
   * \param y2 Составляющая второго вектора положения объекта
   * \param z2 Составляющая второго вектора положения объекта
   * \param vx2 Составляющая второго вектора скорости объекта
   * \param vy2 Составляющая второго вектора скорости объекта
   * \param vz2 Составляющая второго вектора скорости объекта
   * \param t Время интерполяции вектора положения и скорости объекта
   * \param x Составляющая искомого вектора положения объекта
   * \param y Составляющая искомого вектора положения объекта
   * \param z Составляющая искомого вектора положения объекта
   * \param vx Составляющая искомого вектора скорости объекта
   * \param vy Составляющая искомого вектора скорости объекта
   * \param vz Составляющая искомого вектора скорости объекта
   * \return Результат выполнения операции
   */
  static bool interpolatePosVelVector(const T& t1,
                                      const T& x1, const T& y1, const T& z1,
                                      const T& vx1, const T& vy1, const T& vz1,
                                      const T& t2,
                                      const T& x2, const T& y2, const T& z2,
                                      const T& vx2, const T& vy2, const T& vz2,
                                      const T& t,
                                      T& x, T& y, T& z, T& vx, T& vy, T& vz);


/*!
 * \param iearth Модель параметров Земли
 * \param iseason Сезон для расчёта плотности атмосферы
 * \param atm_type Тип атмосферы
 * \param dfront   Величина рубежа
 * \param frx Координата X центра рубежа
 * \param fry Координата Y центра рубежа
 * \param frz Координата Z центра рубежа
 */
  inline CPUBasePropagator(IStaticAtmosphere::ISATypes atm_type, rint32 iseason,
                           const TEarthDatum &iearth, const T &dfront=NAN,
                           const T &frx=NAN, const T &fry=NAN, const T &frz=NAN)
  {
    _datum   = iearth;
    _season  = iseason;
    _atmodel = IStaticAtmosphere::makeInstance(atm_type);
    _deltaFrontier = dfront;
    _frontierX     = frx;
    _frontierY     = fry;
    _frontierZ     = frz;
  }
  virtual inline ~CPUBasePropagator(){ if (_atmodel) delete _atmodel; }


  /*!
   * \brief Расчет переменного баллистического коэффициента
   * \param t1 - время текущего шага
   * \param t2 - время предыдущего шага
   * \param vx1 вектор скорости на текущем шаге
   * \param vy1 вектор скорости на текущем шаге
   * \param vz1 вектор скорости на текущем шаге
   * \param vx2 вектор скорости на предыдущем шаге
   * \param vy2 вектор скорости на предыдущем шаге
   * \param vz2 вектор скорости на предыдущем шаге
   * \param ro плотность атмосферы
   * \return баллистический коэффициент
   */
  T calcVariableGamma(const T &t1, const T &t2,
                      const T &vx1, const T &vy1, const T &vz1,
                      const T &vx2, const T &vy2, const T &vz2, const T &ro);
/*!
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
 * \param gamma Баллистический коэффициент
 * \param gamma1 Баллистический коэффициент производная
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
   * \param gamma Баллистический коэффициент
   * \param gamma1 Баллистический коэффициент производная
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
   * \param atmorange Порог учёта плотности атмосферы
   * \return Результат выполнения операции
   */
  virtual PropagateResult propagate(T &t, const T &dt,
                                    T &x, T &y, T &z,
                                    T &vx, T &vy, T &vz,
                                    T &dvx, T &dvy, T &dvz,
                                    T &pitch, T &raw, T &banking,
                                    T &vpitch, T &vraw, T &vbanking,
                                    T &wx, T &wy, T &wz,
                                    T &dwx, T &dwy, T &dwz,
                                    const T &atmorange=NAN);


/*!
 * \brief Расчёт высоты в базовой системе предсказателя
 * \param x вектор положения в базовой системе предсказателя
 * \param y вектор положения в базовой системе предсказателя
 * \param z вектор положения в базовой системе предсказателя
 */
  virtual T heightInBaseSystem(const T &x, const T &y, const T &z)=0;

  /*!
   * \brief Расчёт параметров точки выведения БР
   * \param Bv Широта точки выведения, [рад]
   * \param Lv Долгота точки выведения, [рад]
   * \param Hv Высота точки выведения, [м]
   * \param Bp Широта точки падения, [рад]
   * \param Lp Долгота точки падения, [рад]
   * \param t Полетное время
   * \param Athimuth Азимут стрельбы
   * \param tetta Угол бросания, [рад]
   * \param x Координата точки выведения в ПГСК, [м]
   * \param y Координата точки выведения в ПГСК, [м]
   * \param z Координата точки выведения в ПГСК, [м]
   * \param vx Элемент вектора скорости в ПГСК, [м\с]
   * \param vy Элемент вектора скорости в ПГСК, [м\с]
   * \param vz Элемент вектора скорости в ПГСК, [м\с]
   * \param iearth Параметры Земли
   * \param opt_angle Флаг расчёта оптимального угла бросания
   * \details Расчет полетного времени, расчет оптимального угла бросания,
   *          расчет составляющих скоростей на начале пассивного участка
   */
  static bool calcInitialVector(const T&Bv, const T&Lv, const T&Hv,
                                const T &Bp, const T &Lp,
                                T &x, T &y, T &z, T &vx, T &vy, T &vz,
                                T &t, T &Athimuth, T &tetta,
                                const TEarthDatum &iearth,
                                bool opt_angle=true);

};


} //namespace RBallistics

} //namespace RMath


#pragma once

#include "math/ballist/cputargeter.h"
#include "math/ballist/cpuballistics.h"

namespace RMath {

using namespace RMatrix;

namespace RBallistics {


/*!
 * \brief Решение задачи Коши для определения начальных условий для изменения
 *        положения объекта с текущего на заданное с учётом максимального
 *        допустимого отклонения.
 */
class CPUTargeterKoshi: public CPUBaseTargeter
{
public:
  inline explicit CPUTargeterKoshi(CPUBasePropagator *iprop):
                                   CPUBaseTargeter(iprop){}
  inline virtual ~CPUTargeterKoshi(){}

  /*!
   * \brief Коррекция требуемых скоростей для попадания в заданную точку с
   *        погрешностью не хуже idelt
   * \param t время
   * \param x Координата текущего положения объекта
   * \param y Координата текущего положения объекта
   * \param z Координата текущего положения объекта
   * \param vx Составляющая скорости объекта в начальном положении
   * \param vy Составляющая скорости объекта в начальном положении
   * \param vz Составляющая скорости объекта в начальном положении
   * \param tx Координата искомого положения объекта
   * \param ty Координата искомого положения объекта
   * \param tz Координата искомого положения объекта
   * \param bk1 БК
   * \param bk2 БК
   * \param idt Шаг интегрирования
   * \param idelt Допустимое отклонение от заданного положения
   * \param maxiter Максимальное количество итераций поиска НУ
   * \param atmorange Порог учёта плотности атмосферы
   * \return Флаг результата операции
   */
  virtual bool adjustConditions(T &t, const T &x, const T &y, const T &z,
                                T &vx, T &vy, T &vz,
                                const T &tx, const T &ty, const T &tz,
                                const T &bk1, const T &bk2,
                                const T& idt, T& idelt, ruint32 maxiter,
                                const T &atmorange=NAN);

};




} //namespace RBallistics

} //namespace RMath


#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif
#include <QtGui>
#include "cpuMSIGMballistics.h"

using namespace RModels::RAtmosphere;

namespace RMath {


namespace RBallistics {


//------------------------------------------------------------------------------
//создание и удаление калькулятора траекторий в МПСК
CPUMPSKPropagator::CPUMPSKPropagator(const T& ib, const T& ih,
                                     IStaticAtmosphere::ISATypes atm_type,
                                     const TEarthDatum &iearth,
                                     const T& ihstab, rint32 iseason,
                                     const T &dfront, const T &frx, const T &fry,
                                     const T &frz):
                                     CPUBasePropagator(atm_type, iseason, iearth,
                                                       dfront, frx, fry, frz)
{
  //расчёт плотности сезонной статической атмосферы на высоте стабилизации
  _rostab = _atmodel->getDensity(ihstab, _season);

  _sinfi = RCore::_sin(ib);
  _cosfi = RCore::_cos(ib);
  T arg  = RCore::_sqrt(1. - _datum.EllipticityEarthDeriv()*_POW2(_sinfi));
#ifdef G_ASSERT_SHEME
   assert(arg!=0.);
#endif
   arg    = _datum.BigSemiAxis()/arg;
  _centerx = 0.;
  _centery = (arg+ih) - _datum.EllipticityEarthDeriv()*arg*_POW2(_sinfi);
  _centerz = -_datum.EllipticityEarthDeriv()*arg*_sinfi*_cosfi;
}

//------------------------------------------------------------------------------
//расчёт высоты объекта над поверхностью эллипсоида по координатам заданным в МПСК
T CPUMPSKPropagator::_heightInMPSK(const T &ycsk, const T &zcsk, const T &dist)
{
  //широта объекта
  T lat = (ycsk*_sinfi + zcsk*_cosfi) / dist;
  T arg = RCore::_sqrt(1. + _datum.Nonsphericity()*_POW2(lat));

  return (RCore::isEqual(arg, 0.)) ?NAN :(dist - _datum.BigSemiAxis()/arg);
}



//---------------------------S252M------------------------------
// Расчет полного ускорения баллистического объекта
// в поле тяготения Земли с учетом атмосферы
//--------------------------------------------------------------
bool CPUMPSKPropagator::accelerations(const T &x, const T &y, const T &z,
                                      const T &vx, const T &vy, const T &vz,
                                      T &dvx, T &dvy, T &dvz,
                                      const T &gamma, const T &gamma1,
                                      const T &atmorange)
{
  //Вектор дальности от объекта к центру Земли в ПГСК
  T  center[3];
  center[0] = _centerx + x;
  center[1] = _centery + y;
  center[2] = _centerz + z;

  //длина вектора из центра Земли на объект
  T distance = RCore::_sqrt(_POW2(center[0])+_POW2(center[1])+_POW2(center[2]));
  if (RCore::isEqual(distance, 0.)) return false;

  //Расчет проекции в-ра дальности объекта из центра Земли  на ось ее вращения
  T fA0 = center[1]*_sinfi + center[2]*_cosfi,
    f1  = fA0/distance;

  //2-й коэф. разложения поля Земли как доля гравитационной постоянной
  static const T CZ = 0.5*_datum.J2()*_POW2(_datum.BigSemiAxis());

  T gravconstR = _datum.GravityConstant()/(_POW3(distance)),

    czR        = CZ/(_POW2(distance)),

    k0         = gravconstR * (czR * (15.*_POW2(f1)-3.) - 1.) +
                 _datum.AngularVelocity2(),

    k1         = fA0 * (6.*gravconstR*czR + _datum.AngularVelocity2());

  //Расчет суммы ускорения гравитации и инерционных ускорений объекта в МПСК
  dvx = k0*center[0] + 2.*_datum.AngularVelocity()*(vy*_cosfi - vz*_sinfi);
  dvy = k0*center[1] - 2.*_datum.AngularVelocity()*vx*_cosfi - k1*_sinfi;
  dvz = k0*center[2] + 2.*_datum.AngularVelocity()*vx*_sinfi - k1*_cosfi;


  //Расчет высоты объекта
  T height = _heightInMPSK(center[1], center[2], distance);
#ifdef G_ASSERT_SHEME
  assert(!RCore::isNan(height));
#else
   if (RCore::isNan(height)) return false;
#endif

  //расчет плотности атмосферы на высоте height
  T density = _atmodel->getDensity(height, _season);
  if ((RCore::isNan(atmorange)==false) && (height>atmorange)) density = 0.;

#ifdef G_ASSERT_SHEME
  assert(!RCore::isNan(density));
#else
   if (RCore::isNan(density)) return false;
#endif

  //проверка необходимости расчета ускорения торможения
  if (density==0.) return true;

  //расчет полного ускорения объекта с учетом влияния атмосферы
  T mvz = RCore::_sqrt(_POW2(vx) + _POW2(vy) + _POW2(vz));

  dvx += ((-0.5)*density*mvz*vx*(gamma + gamma1*RCore::_sqrt(_rostab/density)));
  dvy += ((-0.5)*density*mvz*vy*(gamma + gamma1*RCore::_sqrt(_rostab/density)));
  dvz += ((-0.5)*density*mvz*vz*(gamma + gamma1*RCore::_sqrt(_rostab/density)));

  return true;
}

T CPUMPSKPropagator::heightInBaseSystem(const T &x, const T &y, const T &z)
{
  T X = _centerx + x;
  T Y = _centery + y;
  T Z = _centerz + z;
  T dist = RCore::_sqrt(_POW2(X)+_POW2(Y)+_POW2(Z));
  return _heightInMPSK(Y, Z, dist);
}


//------------------------------------------------------------------------------
//функция расчёта скорости изменения высоты
T CPUMPSKPropagator::altitude_rate(const T &x,  const T &y,  const T &z,
                                   const T &vx, const T &vy, const T &vz)
{
  //Вектор дальности от объекта к центру Земли в ПГСК
  T  center_x = _centerx + x,
     center_y = _centery + y,
     center_z = _centerz + z;

  T rng = RCore::_sqrt(_POW2(center_x) + _POW2(center_y) + _POW2(center_z));
#ifdef G_ASSERT_SHEME
  assert(rng !=0.);
#else
   if (rng==0.) return NAN;
#endif
  rng = 1./rng;
  center_x = center_x*rng;
  center_y = center_y*rng;
  center_z = center_z*rng;

  return (vx*center_x + vy*center_y + vz*center_z);
}

//------------------------------------------------------------------------------
//Прогноз движения объекта в поле Земли на один шаг (в МПСК)
CPUBasePropagator::PropagateResult CPUMPSKPropagator::propagate(T &t, const T &dt,
                                                T &x, T &y, T &z,
                                                T &vx, T &vy, T &vz,
                                                T &dvx, T &dvy, T &dvz,
                                                const T &gamma, const T &gamma1,
                                                const T &atmorange)
{
  if (CPUBasePropagator::propagate(t, dt, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                                   gamma, gamma1, atmorange)!= PR_UNIMPLEMENTED)
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }

  CPUMatrix k0(3, 4),
            k1(3, 4),
            res(10, 1);


  T dt2 = dt/2.;

  //расчет векторных коэффициентов шаг #1
  if (accelerations(x,y,z, vx,vy,vz, dvx,dvy,dvz, gamma, gamma1, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  k0.v(0,0) = vx;
  k0.v(1,0) = vy;
  k0.v(2,0) = vz;
  k1.v(0,0) = dvx;
  k1.v(1,0) = dvy;
  k1.v(2,0) = dvz;
  res[1] = x + k0.v(0,0)*dt2;
  res[2] = y + k0.v(1,0)*dt2;
  res[3] = z + k0.v(2,0)*dt2;
  res[4] = vx + k1.v(0,0)*dt2;
  res[5] = vy + k1.v(1,0)*dt2;
  res[6] = vz + k1.v(2,0)*dt2;

  //шаг #2
  if (accelerations(res[1],res[2],res[3], res[4],res[5],res[6],
                    dvx,dvy,dvz, gamma, gamma1, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  k0.v(0,1) = vx + k1.v(0,0)*dt2;
  k0.v(1,1) = vy + k1.v(1,0)*dt2;
  k0.v(2,1) = vz + k1.v(2,0)*dt2;
  k1.v(0,1) = dvx;
  k1.v(1,1) = dvy;
  k1.v(2,1) = dvz;
  res[1] = x + k0.v(0,1)*dt2;
  res[2] = y + k0.v(1,1)*dt2;
  res[3] = z + k0.v(2,1)*dt2;
  res[4] = vx + k1.v(0,1)*dt2;
  res[5] = vy + k1.v(1,1)*dt2;
  res[6] = vz + k1.v(2,1)*dt2;

  //шаг #3
  if (accelerations(res[1],res[2],res[3], res[4],res[5],res[6],
                    dvx,dvy,dvz, gamma, gamma1, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  k0.v(0,2) = vx + k1.v(0,1)*dt2;
  k0.v(1,2) = vy + k1.v(1,1)*dt2;
  k0.v(2,2) = vz + k1.v(2,1)*dt2;
  k1.v(0,2) = dvx;
  k1.v(1,2) = dvy;
  k1.v(2,2) = dvz;
  res[1] = x + k0.v(0,2)*dt;
  res[2] = y + k0.v(1,2)*dt;
  res[3] = z + k0.v(2,2)*dt;
  res[4] = vx + k1.v(0,2)*dt;
  res[5] = vy + k1.v(1,2)*dt;
  res[6] = vz + k1.v(2,2)*dt;

  //шаг #4
  if (accelerations(res[1],res[2],res[3], res[4],res[5],res[6],
                    dvx,dvy,dvz, gamma, gamma1, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  k0.v(0,3) = vx + k1.v(0,2)*dt;
  k0.v(1,3) = vy + k1.v(1,2)*dt;
  k0.v(2,3) = vz + k1.v(2,2)*dt;
  k1.v(0,3) = dvx;
  k1.v(1,3) = dvy;
  k1.v(2,3) = dvz;
  t += dt;
  x += ( (k0.v(0,0) +2.*k0.v(0,1) + 2.*k0.v(0,2) +  k0.v(0,3))/6.)*dt;
  y += ( (k0.v(1,0) +2.*k0.v(1,1) + 2.*k0.v(1,2) +  k0.v(1,3))/6.)*dt;
  z += ( (k0.v(2,0) +2.*k0.v(2,1) + 2.*k0.v(2,2) +  k0.v(2,3))/6.)*dt;
  vx+= ( (k1.v(0,0) +2.*k1.v(0,1) + 2.*k1.v(0,2) +  k1.v(0,3))/6.)*dt;
  vy+= ( (k1.v(1,0) +2.*k1.v(1,1) + 2.*k1.v(1,2) +  k1.v(1,3))/6.)*dt;
  vz+= ( (k1.v(2,0) +2.*k1.v(2,1) + 2.*k1.v(2,2) +  k1.v(2,3))/6.)*dt;

  if (accelerations(x,y,z, vx,vy,vz, dvx,dvy,dvz, gamma, gamma1, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  return PR_OK;
}


//------------------------------------------------------------------------------
//Прогноз движения объекта в поле Земли до достижения высоты
CPUBasePropagator::PropagateResult
     CPUMPSKPropagator::propagateOnHeight(T &t, const T &dt,
                                          T &x, T &y, T &z,
                                          T &vx, T &vy, T &vz,
                                          T &dvx, T &dvy, T &dvz,
                                          const T &gamma, const T &gamma1,
                                          T &height,
                                          CPUBasePropagator::FlightDirection idir,
                                          const T &maxflighttime, const T &atmorange)
{
  if (CPUBasePropagator::propagate(t, dt, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                                   gamma, gamma1, atmorange)!= PR_UNIMPLEMENTED)
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }

  T accel_vector = NAN;
  const T _criticalAccel = 660.; //крит. значение в-ра ускорения
  const T _criticalAccDiv = 400.;//крит. значение производной в-ра ускорения

  CPUBasePropagator::FlightDirection cur_direction = FD_UNDETERMINED;
  T currentH = heightInBaseSystem(x,y,z),
    oldH     = currentH;
  T at, ax, ay, az, avx, avy, avz;

  if (RCore::isNan(currentH))
  {
    return PR_MATHERROR;
  }

  //полётное время
  T flightTime = 0.;

  bool frontier_cheme = (RCore::isNan(_deltaFrontier)==false) &&
                        (RCore::isNan(_frontierX)==false) &&
                        (RCore::isNan(_frontierY)==false) &&
                        (RCore::isNan(_frontierZ)==false);
  T oldFrontier = NAN;

  while (true)
  {
     at = t; ax = x; ay = y; az = z; avx = vx; avy = vy; avz = vz;

     //прогноз траектории движения объекта на шаг
     CPUBasePropagator::PropagateResult res = propagate(t,dt, x,y,z, vx,vy,vz,
                                                        dvx,dvy,dvz,
                                                        gamma, gamma1, atmorange);
     if (res>CPUBasePropagator::PR_OK_WITH_LIMITATIONS)
     {
       return res;
     }

     //проверка целесообразности прогноза по полётному времени
     flightTime += dt;
     if (flightTime>=maxflighttime)
     {
       return PR_OK_WITH_LIMITATIONS;
     }


     //проверка целесообразности прогноза по критическим параметрам
     T tmp_accel = RCore::_sqrt(_POW2(dvx) + _POW2(dvy) + _POW2(dvz));
     if ( (tmp_accel>_criticalAccel) ||
          ((RCore::isNan(accel_vector)==false) && (accel_vector!=0.) &&
           (_MFASTABS(tmp_accel/accel_vector)>_criticalAccDiv)) )
     {
       t = at; x = ax; y = ay; z = az; vx = avx; vy = avy; vz = avz;
       return PR_OK_WITH_LIMITATIONS;
     }

     //схема работы по рубежу
     if (frontier_cheme)
     {
       T dif  = RCore::_sqrt(_POW2(_frontierX-x) + _POW2(_frontierY-y) +
                             _POW2(_frontierZ-z));

       oldFrontier = dif;
       if ( (RCore::isNan(oldFrontier)==false) && (RCore::isNan(dif)==false) &&
            (dif > oldFrontier) &&  (dif > _deltaFrontier))
       {
         return PR_FRONTIER;
       }
     }

     accel_vector = tmp_accel;

     //определение высоты и направления траектории
     currentH = heightInBaseSystem(x,y,z);
     if (RCore::isNan(currentH))
     {
       return PR_MATHERROR;
     }
     cur_direction = (oldH>currentH)  ?FD_DESCENDING
                     :(oldH<currentH) ?FD_ASCENDING :FD_UNDETERMINED;
     oldH = currentH;

     //уточнение попадания в пороговое значение высоты
     if ( ((cur_direction==idir) || (idir==FD_UNDETERMINED)) &&
          //на взлёте превысили порог высоты
          ( ((cur_direction==FD_ASCENDING) && (currentH>=height)) ||
          //при снижении стали меньше нижнего порога
            ((cur_direction==FD_DESCENDING) && (currentH<=height)) ) )
     {

        T smallstep = -0.00001;
        while ( ((cur_direction==FD_ASCENDING) && (currentH>height)) ||
                ((cur_direction==FD_DESCENDING) && (currentH<height)) )
        {
            res = propagate(t,smallstep, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                            gamma, gamma1, atmorange);
            if (res!=PR_OK)
            {
              return PR_MATHERROR;
            }
            currentH = this->heightInBaseSystem(x,y,z);
        }
        return res;
     }
  }

  return PR_OK;
}





} //namespace RBallistics

} //namespace RMath

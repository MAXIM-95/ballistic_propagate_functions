#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cpu45ballistics.h"

namespace RMath {



namespace RBallistics {


//------------------------------------------------------------------------------
//Коэффициент №1                                                               |
T CPU45Propagator::K1(const T &x, const T &y, const T &z,
                      const T &Vx, const T &Vy, const T &Vz,
                      const T &dt, const T &bk, const T &atmorange)
{
  T   Wz = _datum.AngularVelocity(),
      r  =  RCore::_sqrt(_POW2(x)  + _POW2(y)  + _POW2(z)),
      rV =  RCore::_sqrt(_POW2(Vx) + _POW2(Vy) + _POW2(Vz));

#ifdef G_ASSERT_SHEME
  assert(r!=0.);
#else
  if (r==(T)0.) return NAN;
#endif

  T   b   = TEarthDatum::meanRadius()/_POW3(r),
      c   = 1.5*_datum.a20()*_POW2(TEarthDatum::meanRadius())/_POW2(r);
  T   d   = z/r;
      d   = 5.*_POW2(d);
  T   a   = b*(_datum.a00() + c*(d - 1.)),
      hgt = heightInPGSK(r, z, _datum);
  //расчет плотности атмосферы на высоте height
  T density = _atmodel->getDensity(hgt, _season);
  if ((RCore::isNan(atmorange)==false) && (hgt>atmorange)) density = 0.;

  return ((_POW2(Wz) - a)*x + 2.*Wz*Vy - 0.5*bk*density*Vx*rV)*dt;
}


//------------------------------------------------------------------------------
//Коэффициент №2                                                               |
T CPU45Propagator::K2(const T &x, const T &y, const T &z,
                      const T &Vx, const T &Vy, const T &Vz,
                      const T &dt, const T &bk, const T &atmorange)
{
  T   Wz = _datum.AngularVelocity(),
      r  =  RCore::_sqrt(_POW2(x)  + _POW2(y)  + _POW2(z)),
      rV =  RCore::_sqrt(_POW2(Vx) + _POW2(Vy) + _POW2(Vz));

#ifdef G_ASSERT_SHEME
  assert(r!=0.);
#else
  if (r==(T)0.) return NAN;
#endif

  T   b   = TEarthDatum::meanRadius()/_POW3(r),
      c   = 1.5*_datum.a20()*_POW2(TEarthDatum::meanRadius())/_POW2(r);
  T   d   = z/r;
      d   = 5.*_POW2(d);
  T   a   = b*(_datum.a00() + c*(d - 1.)),
      hgt = heightInPGSK(r, z, _datum);
  //расчет плотности атмосферы на высоте height
  T density = _atmodel->getDensity(hgt, _season);
  if ((RCore::isNan(atmorange)==false) && (hgt>atmorange)) density = 0.;

  return((_POW2(Wz) - a)*y - 2*Wz*Vx - 0.5*bk*density*Vy*rV)*dt;
}

//------------------------------------------------------------------------------
//Коэффициент №3                                                               |
T CPU45Propagator::K3(const T &x, const T &y, const T &z,
                      const T &Vx, const T &Vy, const T &Vz,
                      const T &dt, const T &bk, const T &atmorange)
{
  T   r  =  RCore::_sqrt(_POW2(x)  + _POW2(y)  + _POW2(z)),
      rV =  RCore::_sqrt(_POW2(Vx) + _POW2(Vy) + _POW2(Vz));

#ifdef G_ASSERT_SHEME
  assert(r!=0.);
#else
  if (r==(T)0.) return NAN;
#endif

  T   b   = TEarthDatum::meanRadius()/_POW3(r),
      c   = 1.5*_datum.a20()*_POW2(TEarthDatum::meanRadius())/_POW2(r);
  T   d   = z/r;
      d   = 5.*_POW2(d);
  T   a   = b*(_datum.a00() + c*(d - 1.)),
      hgt = heightInPGSK(r, z, _datum);
  //расчет плотности атмосферы на высоте height
  T density = _atmodel->getDensity(hgt, _season);
  if ((RCore::isNan(atmorange)==false) && (hgt>atmorange)) density = 0.;

  return((2.*b*c - a)*z - 0.5*bk*density*Vz*rV)*dt;
}


//------------------------------------------------------------------------------
//высота над поверхностью эллипсоида
T CPU45Propagator::heightInBaseSystem(const T &x, const T &y, const T &z)
{
  return heightInPGSK(x, y, z, _datum);
}



//------------------------------------------------------------------------------
//создание и удаление калькулятора траекторий в ПГСК
CPU45Propagator::CPU45Propagator(IStaticAtmosphere::ISATypes atm_type,
                                 rint32 iseason, const TEarthDatum &iearth,
                                 const T &dfront, const T &frx, const T &fry,
                                 const T &frz):
                                 CPUBasePropagator(atm_type, iseason, iearth,
                                                   dfront, frx, fry, frz)
{

}


//------------------------------------------------------------------------------
//функция расчёта скорости изменения высоты
T CPU45Propagator::altitude_rate(const T &x,  const T &y,  const T &z,
                                   const T &vx, const T &vy, const T &vz)
{
  //Вектор дальности от объекта к центру Земли в ПГСК
  T  center_x = x,
     center_y = y,
     center_z = z;

  T rng = RCore::_sqrt(_POW2(center_x) + _POW2(center_y) + _POW2(center_z));
#ifdef G_ASSERT_SHEME
  assert(rng !=0.);
#else
   if (rng==0.) return NAN;
#endif
  rng = 1./rng;
  center_x *= rng;
  center_y *= rng;
  center_z *= rng;

  return (vx*center_x + vy*center_y + vz*center_z);
}

//------------------------------------------------------------------------------
//Прогноз движения объекта в поле Земли на один шаг
CPUBasePropagator::PropagateResult
                 CPU45Propagator::propagate(T &t, const T &dt, T &x, T &y, T &z,
                                T &vx, T &vy, T &vz, T &dvx, T &dvy, T &dvz,
                                const T &gamma, const T &gamma1, const T &atmorange)
{
  if (CPUBasePropagator::propagate(t, dt, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                                   gamma, gamma1, atmorange)!= PR_UNIMPLEMENTED)
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }

  dvx = dvy = dvz = NAN;

  T K11 = K1(x, y, z, vx, vy, vz, dt, gamma, atmorange),
    K12 = K2(x, y, z, vx, vy, vz, dt, gamma, atmorange),
    K13 = K3(x, y, z, vx, vy, vz, dt, gamma, atmorange),
    K14 = vx*dt,
    K15 = vy*dt,
    K16 = vz*dt;
  if (RCore::isNan(K11) || RCore::isNan(K12) || RCore::isNan(K13))
  {
    return PR_MATHERROR;
  }

  T kx  = 0.5*K14,
    ky  = 0.5*K15,
    kz  = 0.5*K16,
    kvx = 0.5*K11,
    kvy = 0.5*K12,
    kvz = 0.5*K13;

  T K21 = K1(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K22 = K2(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K23 = K3(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K24 = (vx + kvx)*dt,
    K25 = (vy + kvy)*dt,
    K26 = (vz + kvz)*dt;
  if (RCore::isNan(K21) || RCore::isNan(K22) || RCore::isNan(K23))
  {
    return PR_MATHERROR;
  }

  kx  = 0.5*K24;
  ky  = 0.5*K25;
  kz  = 0.5*K26;
  kvx = 0.5*K21;
  kvy = 0.5*K22;
  kvz = 0.5*K23;

  T K31 = K1(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K32 = K2(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K33 = K3(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K34 = (vx + kvx)*dt,
    K35 = (vy + kvy)*dt,
    K36 = (vz + kvz)*dt;
  if (RCore::isNan(K31) || RCore::isNan(K32) || RCore::isNan(K33))
  {
    return PR_MATHERROR;
  }

  kx  = K34;
  ky  = K35;
  kz  = K36;
  kvx = K31;
  kvy = K32;
  kvz = K33;

  T K41 = K1(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K42 = K2(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K43 = K3(x + kx, y + ky, z + kz, vx + kvx, vy + kvy, vz + kvz, dt, gamma, atmorange),
    K44 = (vx + kvx)*dt,
    K45 = (vy + kvy)*dt,
    K46 = (vz + kvz)*dt;
  if (RCore::isNan(K41) || RCore::isNan(K42) || RCore::isNan(K43))
  {
    return PR_MATHERROR;
  }

  t  += dt;
  x  += (K14 + 2.*K24 + 2.*K34 + K44)/6.;
  y  += (K15 + 2.*K25 + 2.*K35 + K45)/6.;
  z  += (K16 + 2.*K26 + 2.*K36 + K46)/6.;
  vx += (K11 + 2.*K21 + 2.*K31 + K41)/6.;
  vy += (K12 + 2.*K22 + 2.*K32 + K42)/6.;
  vz += (K13 + 2.*K23 + 2.*K33 + K43)/6.;

  return PR_OK;
}




//------------------------------------------------------------------------------
//Прогноз движения объекта в поле Земли до достижения высоты
CPUBasePropagator::PropagateResult
     CPU45Propagator::propagateOnHeight(T &t, const T &dt,
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

     //схема работы по рубежу
     if (frontier_cheme)
     {
       T dif  = RCore::_sqrt(_POW2(_frontierX-x) + _POW2(_frontierY-y) +
                             _POW2(_frontierZ-z));
       if (RCore::isNan(oldFrontier)) oldFrontier = dif;
       else if ( (dif > oldFrontier) && (dif > _deltaFrontier))
       {
         return PR_FRONTIER;
       }
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
            currentH = heightInBaseSystem(x,y,z);
        }
        return res;
     }
  }

  return PR_OK;
}



} //namespace RBallistics

} //namespace RMath

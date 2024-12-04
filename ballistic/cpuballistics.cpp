#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cpuballistics.h"

using namespace RModels::RAtmosphere;

namespace RMath {


namespace RBallistics {




CPUBasePropagator::PropagateResult
                        CPUBasePropagator::propagate(T &t, const T &dt,
                                                     T &x, T &y, T &z,
                                                     T &vx, T &vy, T &vz,
                                                     T &dvx, T &dvy, T &dvz,
                                                     const T &gamma,
                                                     const T &gamma1,
                                                     const T &)
{
    if (RCore::isNan(t)     || RCore::isNan(dt)   ||
        RCore::isNan(x)     || RCore::isNan(y)    || RCore::isNan(z)  ||
        RCore::isNan(vx)    || RCore::isNan(vy)   || RCore::isNan(vz) ||
        RCore::isNan(dvx)   || RCore::isNan(dvy)  || RCore::isNan(dvz)||
        RCore::isNan(gamma) || RCore::isNan(gamma1))
    {
      return PR_INVALID_INPUT_PARAMETERS;
    }

    return PR_UNIMPLEMENTED;
}


CPUBasePropagator::PropagateResult
                   CPUBasePropagator::propagate(T &t, const T &dt,
                                                T &x, T &y, T &z,
                                                T &vx, T &vy, T &vz,
                                                T &dvx, T &dvy, T &dvz,
                                                T &pitch, T &raw, T &banking,
                                                T &vpitch, T &vraw, T &vbanking,
                                                T &wx, T &wy, T &wz,
                                                T &dwx, T &dwy, T &dwz,
                                                const T &)
{
    if (RCore::isNan(t)     || RCore::isNan(dt)   ||
        RCore::isNan(x)     || RCore::isNan(y)    || RCore::isNan(z)  ||
        RCore::isNan(vx)    || RCore::isNan(vy)   || RCore::isNan(vz) ||
        RCore::isNan(dvx)   || RCore::isNan(dvy)  || RCore::isNan(dvz)||
        RCore::isNan(pitch) || RCore::isNan(raw)  || RCore::isNan(banking) ||
        RCore::isNan(vpitch)|| RCore::isNan(vraw) || RCore::isNan(vbanking)||
        RCore::isNan(wx )   || RCore::isNan(wy)   || RCore::isNan(wz)||
        RCore::isNan(dwx)   || RCore::isNan(dwy)  || RCore::isNan(dwz) )
    {
      return PR_INVALID_INPUT_PARAMETERS;
    }

    return PR_UNIMPLEMENTED;
}


CPUBasePropagator::PropagateResult
        CPUBasePropagator::propagateOnHeight(T &t, const T &dt,
                                             T &x, T &y, T &z,
                                             T &vx, T &vy, T &vz,
                                             T &dvx, T &dvy, T &dvz,
                                             const T &gamma, const T &gamma1,
                                             T &height, FlightDirection,
                                             const T &maxflighttime,
                                             const T &)
{
  if (RCore::isNan(t)     || RCore::isNan(dt)   ||
      RCore::isNan(x)     || RCore::isNan(y)    || RCore::isNan(z)  ||
      RCore::isNan(vx)    || RCore::isNan(vy)   || RCore::isNan(vz) ||
      RCore::isNan(dvx)   || RCore::isNan(dvy)  || RCore::isNan(dvz)||
      RCore::isNan(gamma) || RCore::isNan(gamma1)  || RCore::isNan(height) ||
      RCore::isNan(maxflighttime) )
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }

  return PR_UNIMPLEMENTED;
}



//------------------------------------------------------------------------------
//Расчёт высоты объекта над эллипсоидом
T CPUBasePropagator::heightInPGSK(const T &r, const T &z,
                                  const TEarthDatum &datum)
{
#ifdef G_ASSERT_SHEME
  assert(r!=0.);
#else
  if (r==0.) return NAN;
#endif
  const T _e2 = _POW2(datum.Eccentricity());
  T H = 1.0 - _e2*_POW2(RCore::_cos(RCore::_asin(z/r)));
#ifdef G_ASSERT_SHEME
  assert(H!=0.);
#else
  if (H==0.) return NAN;
#endif

  return r - datum.BigSemiAxis()*RCore::_sqrt((1.0 - _e2)/H);
}



const char *CPUBasePropagator::PropagateResultDescription(PropagateResult in)
{
    switch (in)
    {
      case PR_OK:                       return "Штатное завершение";
      case PR_OK_WITH_LIMITATIONS:      return "Прогноз прерван по критическому условию";
      case PR_INVALID_INPUT_PARAMETERS: return "Некорректные входные параметры";
      case PR_MATHERROR:                return "Прерван - ошибка в вычислениях";
      case PR_FRONTIER:                 return "Прерван - выход за указанный рубеж";

      default:                          return "Отсутствует реализация";
    }
}

T CPUBasePropagator::heightInPGSK(const T &x, const T &y, const T &z,
                                  const TEarthDatum &datum)
{
  T r = RCore::_sqrt(_POW2(x) + _POW2(y) + _POW2(z));

  return heightInPGSK(r, z, datum);
}

//------------------------------------------------------------------------------
//Интерполяция положения объекта
bool CPUBasePropagator::interpolatePosVelVector(const T &t1,
                                                const T &x1, const T &y1,
                                                const T &z1, const T &vx1,
                                                const T &vy1, const T &vz1,
                                                const T &t2,
                                                const T &x2, const T &y2,
                                                const T &z2, const T &vx2,
                                                const T &vy2, const T &vz2,
                                                const T &t, T &x, T &y, T &z,
                                                T &vx, T &vy, T &vz)
{
  if (t==t1)
  {
    x = x1; y = y1; z = z1; vx = vx1; vy = vy1; vz = vz1;
    return true;
  }
  else if (t==t2)
  {
    x = x2; y = y2; z = z2; vx = vx2; vy = vy2; vz = vz2;
    return true;
  }

  T tau   = t2 - t1,
    taui  = t - t1;

  if (RCore::isEqual(tau, 0.))  return false;

  T tau2  = _POW2(tau),
    tau3  = _POW3(tau),
    taui2 = _POW2(taui),
    taui3 = _POW3(taui);
  T dx = x2 - x1,
    dy = y2 - y1,
    dz = z2 - z1,
   k2x = (3.*dx - tau*(vx2 + 2.*vx1))/tau2,
   k2y = (3.*dy - tau*(vy2 + 2.*vy1))/tau2,
   k2z = (3.*dz - tau*(vz2 + 2.*vz1))/tau2,
   k3x = (-2.*dx + tau*(vx2 + vx1))/tau3,
   k3y = (-2.*dy + tau*(vy2 + vy1))/tau3,
   k3z = (-2.*dz + tau*(vz2 + vz1))/tau3;

   x  = x1 + vx1*taui + k2x*taui2 + k3x*taui3;
   y  = y1 + vy1*taui + k2y*taui2 + k3y*taui3;
   z  = z1 + vz1*taui + k2z*taui2 + k3z*taui3;
   vx = vx1 + 2.*k2x*taui + 3.*k3x*taui2;
   vy = vy1 + 2.*k2y*taui + 3.*k3y*taui2;
   vz = vz1 + 2.*k2z*taui + 3.*k3z*taui2;

   return true;
}


//-------------------------------------------------------------------------------
//расчёт переменного БК
T CPUBasePropagator::calcVariableGamma(const T &t1, const T &t2, const T &vx1,
                                       const T &vy1, const T &vz1, const T &vx2,
                                       const T &vy2, const T &vz2, const T &ro)
{
  T V1 = RCore::_sqrt(_POW2(vx1)+_POW2(vy1)+_POW2(vz1));
  T V2 = RCore::_sqrt(_POW2(vx2)+_POW2(vy2)+_POW2(vz2));
  T q = t1-t2;
#ifdef G_ASSERT_SHEME
   assert(q!=(T)0.);
#else
   if (q==(T)0.) return NAN;
#endif
  T a = (V1 - V2)/q;
  q = 0.5*ro*V2;
#ifdef G_ASSERT_SHEME
   assert(q!=(T)0.);
#else
   if (q==(T)0.) return NAN;
#endif
  return a/q;
}





//------------------------------------------------------------------------------
// расчёт параметров точки выведения БР
bool CPUBasePropagator::calcInitialVector(const T&Bv, const T&Lv, const T&Hv,
                                          const T &Bp, const T &Lp,
                                          T &x, T &y, T &z, T &vx, T &vy, T &vz,
                                          T &t, T &Athimuth,  T &tetta,
                                          const TEarthDatum &iearth,
                                          bool opt_angle)
{
   const T b0 = iearth.GravityConstant();
   t     = 700.0;
   T  angv  = (opt_angle) ?_MDEG2RAD(45.) :tetta,
      r0    = TEarthDatum::meanRadius() + Hv,
      r0_   = r0/TEarthDatum::meanRadius(),
      V0    = 0.,
      p_t   = t,
      Fp,
      sinbv = RCore::_sin(Bv),
      cosbv = RCore::_cos(Bv),
      sinbp = RCore::_sin(Bp),
      cosbp = RCore::_cos(Bp);

#ifdef G_ASSERT_SHEME
   assert(r0_!=(T)0.);
#else
   if (r0_==(T)0.) return false;
#endif

   //рассчитываем полётное время
   do
   {
     T lamd = Lp + iearth.AngularVelocity()*t;
       Fp   = RCore::_acos(sinbv*sinbp + RCore::_cos(lamd - Lv)*cosbv*cosbp);
     T tfp2 = RCore::_tan(0.5*Fp),
       tang = RCore::_tan(angv),
       cosa = RCore::_cos(angv),
       nu   =  (r0_ + 1.)*_POW2(tfp2) + 2.*tfp2*tang + r0_ - 1.;

#ifdef G_ASSERT_SHEME
     assert(nu!=(T)0.);
#else
     if (nu==(T)0.) return false;
#endif
       nu   = (2.*(1. + _POW2(tang))*_POW2(tfp2))/nu;
     T p    = r0*nu*_POW2(cosa),
       e_   = 1. - (2. - nu)*nu*_POW2(cosa);

#ifdef G_ASSERT_SHEME
     assert(e_>=(T)0.);
#else
     if (e_<(T)0.) return false;
#endif
     e_     = sqrt(e_);

     T nup  = 2. + (nu - 2.)/r0_;
     p_t  = t;

#ifdef G_ASSERT_SHEME
     assert(nup!=(T)0.);
#else
     if (nup==(T)0.) return false;
#endif
     T tettap    = nu/nup;
#ifdef G_ASSERT_SHEME
     assert(tettap>=0.);
#else
     if (tettap<(T)0.) return false;
#endif
     tettap    = sqrt(tettap);
     tettap    = r0_*cosa*tettap;
#ifdef G_ASSERT_SHEME
     assert((tettap>=-1.) && (tettap<=1.));
#else
     if ( (tettap<(T)-1.) || (tettap>(T)1.) ) return false;
#endif
     tettap    = -RCore::_acos(tettap);

#ifdef G_ASSERT_SHEME
     assert(e_ >= 0.);
     assert(p  >= 0.);
#else
     if (e_<(T)0.) return false;
     if (p<(T)0.)  return false;
#endif

     t = 1. - _POW2(e_);
#ifdef G_ASSERT_SHEME
     assert(t>=0.);
#else
     if (t<(T)0.) return false;
#endif
     t = p/t;
     t *= RCore::_sqrt(p/b0);

     T arg1 = (1. - nup)/e_;
     T arg2 = (1. - nu)/e_;
#ifdef G_ASSERT_SHEME
     assert((arg1>=-1.) && (arg1<=1.));
     assert((arg2>=-1.) && (arg2<=1.));
#else
     if ( (arg1<(T)-1.) || (arg1>(T)1.) ) return false;
     if ( (arg2<(T)-1.) || (arg2>(T)1.) ) return false;
#endif
     arg1 = RCore::_acos(arg1);
     arg2 = RCore::_acos(arg2);
     T arg = 1. - _POW2(e_);
#ifdef G_ASSERT_SHEME
     assert(arg>=0.);
#else
     if (arg<(T)0.) return false;
#endif
     arg = RCore::_sqrt(arg);
#ifdef G_ASSERT_SHEME
     assert(arg>=0.);
#else
     if (arg<(T)0.) return false;
#endif
     t *= (arg1 + arg2)/arg + tang - RCore::_tan(tettap);

     V0 = nu*b0/r0;
#ifdef G_ASSERT_SHEME
     assert(V0>=0.);
#else
     if (V0<(T)0.) return false;
#endif
     V0 = RCore::_sqrt(V0);
     Athimuth = asin(RCore::_sin(lamd - Lv)*cosbp/RCore::_sin(Fp));
  }
  while(_MFASTABS(t - p_t)>=10.);

  //расчёт оптимального угла бросания
  if (opt_angle) tetta = _MDEG2RAD(45.0) - Fp/4.0;

  //рассчитываем составляющие скоростей на начало пассивного участка
  T cosangv = RCore::_cos(angv);
  CPUMatrix MV0(1, 3);
  MV0[0] = V0*cosangv*RCore::_cos(Athimuth);
  MV0[1] = V0*RCore::_sin(angv);
  MV0[2] = V0*cosangv*RCore::_sin(Athimuth);

  T sinlv = RCore::_sin(Lv),
    coslv = RCore::_cos(Lv);
  CPUMatrix MC(3,3);
  MC.v(0,0) = -sinbv*coslv;
  MC.v(1,0) = -sinbv*sinlv;
  MC.v(2,0) = cosbv;
  MC.v(0,1) = cosbv*coslv;
  MC.v(1,1) = cosbv*sinlv;
  MC.v(2,1) = sinbv;
  MC.v(0,2) = -sinlv;
  MC.v(1,2) = coslv;
  MC.transpose(false);
  CPUMatrix MV = MC*MV0;

  if (!geodez2pgsk(Bv, Lv, Hv, x, y, z, iearth))
  {
    return false;
  }
  vx = MV[0] + iearth.AngularVelocity()*y;
  vy = MV[1] - iearth.AngularVelocity()*x;
  vz = MV[2];


  return true;
}


} //namespace RBallistics

} //namespace RMath

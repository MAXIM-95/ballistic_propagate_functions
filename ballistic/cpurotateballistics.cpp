#include <string.h>
#include <stdlib.h>
#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cpurotateballistics.h"

using namespace RModels::RAtmosphere;

double _uniform()
{
   // srand(time(NULL));
    int r = rand();
    return (double)r / RAND_MAX;
}

double _uniformMINtoMAX(double min, double max)
{
    return _uniform() * (max - min) + min;
}

double _exponential(double lambda)
{
    return -log(_uniform()) / lambda;
}

double _normalhelper()
{
    double c = sqrt(2 * M_E / M_PI);
    double t = _exponential(1.0);

    while(t > (sqrt(2 / M_PI) * exp(t - t * t / 2) / c))
    {
        t = _exponential(1.0);
    }

    if(rand() % 2 == 0)
    {
        t = -1 * t;
    }

    return 2 * t;
}

double _normal(double mu, double sigma)
{
    return sigma * _normalhelper() + mu;
}

namespace RMath {


namespace RBallistics {

CPUElemParams::CPUElemParams(const T &iwx, const T &iwy, const T &iwz)
{
  memset((void *)this, 0, sizeof(struct CPUElemParams));
  wx = iwx;
  wy = iwy;
  wz = iwz;
  //==================
 // wx = _normal(0.,0.2);
 // wy = _normal(0.,0.2);
 // wz = _normal(18.78,0.);
}

CPUElemParams &CPUElemParams::operator=(const CPUElemParams &other)
{
  if (this==&other) return *this;

  memcpy((void *)this, &other, sizeof(struct CPUElemParams));
  return *this;
}


//------------------------------------------------------------------------------
//Расстояние от носика до центра тяжести вдоль продольной оси
T CPUElemParams::centerOfGravityShift()
{
  T s1 = 0.5*M_PI*_POW2(nose_radius),
    c1 = nose_radius*(1. - 4./(3.*M_PI)),
    h  = length - nose_radius,
    s2 = h*(nose_radius*0.5*base_diameter),
    c2 = base_diameter + 2.*nose_radius;
  //расчёт угловых ускорений
#ifdef G_ASSERT_SHEME
  assert(c2!=0.);
#else
  if (c2==0.) return NAN;
#endif
  c2 = nose_radius + h/3.*(base_diameter/c2 + 1.);
  T res = s1 + s2;
#ifdef G_ASSERT_SHEME
  assert(res!=0.);
#else
  if (res==0.) return NAN;
#endif

  return (s1*c1 + s2*c2)/res;
}


//------------------------------------------------------------------------------
//Расчёт начальных углов на входе в атмосферу
bool CPUElemParams::calculateInitialAngles(const T &vx, const T &vy,
                                           const T &vz,
                                           T &pitch, T &raw, T &banking,
                                           T &owx, T &owy, T &owz)
{
  T v = RCore::_sqrt(_POW2(vx) + _POW2(vy) + _POW2(vz));
  T vxz = RCore::_sqrt(_POW2(vx) + _POW2(vz));

  owx = wx;
  owy = wy;
  owz = wz;

#ifdef G_ASSERT_SHEME
  assert(v!=0.);
#else
  if (v==0.) return false;
#endif

  pitch = _MSIGN(vy)*RCore::_acos(vxz/v) + attack_angle;
  T tmp = RCore::isEqual(vxz, 0.) ?0. :vx/vxz;
  raw = (RCore::isEqual(vxz, 0.)) ?0.
        :(vz<=0.) ?0.5*M_PI - RCore::_asin(tmp) :1.5*M_PI + RCore::_asin(tmp);
  banking = 0.;

  return true;
}








//------------------------------------------------------------------------------
//создание и удаление калькулятора траекторий в МПСК (Моисеевский)
CPUMPSKRotatePropagator::CPUMPSKRotatePropagator(const T& ib, const T& ih,
                                                 const T &imaxheight,
                                     IStaticAtmosphere::ISATypes atm_type,
                                     const TEarthDatum &iearth,
                                     rint32 iseason, T icurr,
                                     const T &dfront, const T &frx,
                                     const T &fry, const T &frz):
                                     CPUBasePropagator(atm_type, iseason, iearth,
                                                       dfront, frx, fry, frz)
{
  _currency  = icurr;
  _maxheight = imaxheight;
  _params    = nullptr;
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
T CPUMPSKRotatePropagator::_heightInMPSK(const T &ycsk, const T &zcsk, const T &dist)
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
bool CPUMPSKRotatePropagator::accelerations(const T &x, const T &y, const T &z,
                                            const T &vx, const T &vy, const T &vz,
                                            T &dvx, T &dvy, T &dvz,
                                            const T &pitch, const T &raw, const T &banking,
                                            T &vpitch, T &vraw, T &vbanking,
                                            T &wx, T &wy, T &wz,
                                            T &dwx, T &dwy, T &dwz, bool inatm,
                                            T &hght, const T &atmorange)
{
  //Вектор дальности от объекта к центру Земли в ПГСК
  T  center[3];
  center[0] = _centerx + x;
  center[1] = _centery + y;
  center[2] = _centerz + z;

  //длина вектора из центра Земли на объект
  T distance = RCore::_sqrt(_POW2(center[0])+_POW2(center[1])+_POW2(center[2]));
  if (RCore::isEqual(distance, 0.)) return false;

  //Расчет проекции в-ра дальности объекта из центра Земли на ось ее вращения
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

  hght = NAN;

  if (inatm==false) return true;

  //Расчет высоты объекта
  hght = _heightInMPSK(center[1], center[2], distance);
#ifdef G_ASSERT_SHEME
  assert(RCore::isNan(hght)==false);
#else
   if (RCore::isNan(hght)) return false;
#endif

  //расчет плотности атмосферы на высоте height
  T density = _atmodel->getDensity(hght, _season);
  if ((RCore::isNan(atmorange)==false) && (hght>atmorange)) density = 0.;

#ifdef G_ASSERT_SHEME
  assert(!RCore::isNan(density));
#else
   if (RCore::isNan(density)) return false;
#endif

  //проверка необходимости расчета ускорения торможения
  if (density==0.) return true;

  //расчет полного ускорения объекта с учетом влияния атмосферы
  T mvz = RCore::_sqrt(_POW2(vx) + _POW2(vy) + _POW2(vz));
#ifdef G_ASSERT_SHEME
  assert(mvz !=0.);
#else
   if (mvz==0.) return false;
#endif
  T qt = 0.5*density*_params->midel*_POW2(mvz);

  CPUMatrix *ssk = MPSK2LSK(pitch, raw, banking);
  if (ssk==nullptr) return false;
  CPUMatrix V(1, 3);
  V[0] = vx;
  V[1] = vy;
  V[2] = vz;

  //проверка на восходящее движение в атмосфере
  CPUMatrix Mcenter(3, 1);
  Mcenter[0] = center[0];
  Mcenter[1] = center[1];
  Mcenter[2] = center[2];

  V = (*ssk)*V;

  const T delt = 1e-12;
  T vyz = RCore::_sqrt(_POW2(V[1]) + _POW2(V[2]));
  T attAng = (vyz>delt) ?RCore::_acos(V[0]/mvz) :0.;
  T ang1   = (vyz>delt) ?V[2]/vyz               :0.;
  T ang2   = (vyz>delt) ?V[1]/vyz               :1.;
  T Ctau   = 0.03616667 + 0.29671029*attAng - 0.01406917*_POW2(attAng);
  T Cn     = 2.03972975*attAng;
  T Fn     = -qt*Cn;

#ifdef G_ASSERT_SHEME
  assert(_params->weight !=0.);
#else
   if (_params->weight==0.) return false;
#endif
  ssk->transpose(false);
  T mass = (1./_params->weight);

  CPUMatrix Fxyz(1, 3);
  Fxyz[0]  = -qt*Ctau;
  Fxyz[1]  = Fn*ang2;
  Fxyz[2]  = Fn*ang1;
  CPUMatrix Fmass(Fxyz);
  Fmass[0] = Fmass[0]*mass;
  Fmass[1] = Fmass[1]*mass;
  Fmass[2] = Fmass[2]*mass;
  Fmass = *ssk * Fmass;

  dvx += Fmass[0];
  dvy += Fmass[1];
  dvz += Fmass[2];
  delete ssk;

  T sa = RCore::_sin(banking),
  ca = RCore::_cos(banking);
  vpitch = wy*sa + wz*ca;
  vraw   = RCore::_cos(pitch);
#ifdef G_ASSERT_SHEME
  assert(vraw!=0.);
#else
   if (vraw==0.) return false;
#endif
  vraw = 1./vraw*(wy*ca - wz*sa);
  vbanking = wx - sa*vraw;

  //расчёт восстанавливающего момента
  CPUMatrix Ma(1,3);
  if (attAng>0.)
  {
    static T dx = 0.075*_params->length;
    Ma[1] =  Fxyz[2]*dx;
    Ma[2] = -Fxyz[1]*dx;
  }

  //расчёт демпфирующего момента
  qt = 0.5*density*_params->midel*mvz;
  CPUMatrix Cxyz(1, 3);
  Cxyz[0] = 0.1;
  Cxyz[1] = 0.17;
  Cxyz[2] = 0.17;

  const T Myz = -Cxyz[1]*qt*_POW2(_params->length);
  CPUMatrix Md(1,3);
  Md[0] = -Cxyz[0]*qt*_POW2(_params->base_diameter)*wx;
  Md[1] = Myz*wy;
  Md[2] = Myz*wz;

  CPUMatrix Mad = Ma + Md;

  //расчёт угловых ускорений
#ifdef G_ASSERT_SHEME
  assert(_params->inert_moment_x!=0.);
  assert(_params->inert_moment_y!=0.);
  assert(_params->inert_moment_z!=0.);
#else
  if (_params->inert_moment_x==0.) return false;
  if (_params->inert_moment_y==0.) return false;
  if (_params->inert_moment_z==0.) return false;
#endif

  dwx = (Mad[0] + (_params->inert_moment_y - _params->inert_moment_z)*wy*wz )/
                   _params->inert_moment_x;
  dwy = (Mad[1] + (_params->inert_moment_z - _params->inert_moment_x)*wx*wz )/
                   _params->inert_moment_y;
  dwz = (Mad[2] + (_params->inert_moment_x - _params->inert_moment_y)*wx*wy )/
                   _params->inert_moment_z;

  return true;
}




//------------------------------------------------------------------------------
//Преобразование МПСК ==> Связанная СК
CPUMatrix *CPUMPSKRotatePropagator::MPSK2LSK(const T &pitch, const T &raw,
                                             const T &banking)
{
  const T c0 = RCore::_cos(pitch),
          c1 = RCore::_cos(raw),
          c2 = RCore::_cos(banking),
          s0 = RCore::_sin(pitch),
          s1 = RCore::_sin(raw),
          s2 = RCore::_sin(banking);
  T arr[] =
  {
      c0*c1,               s0,             -s1*c0,
     -c1*s0*c2 + s1*s2,    c0*c2,           s0*s1*c2 + c1*s2,
      c1*s0*s2 + s1*c2,    -c0*s2,          -s0*s1*s2 + c1*c2
  };

  return CPUMatrix::from2DArray(arr, 3, 3, true);
}



//------------------------------------------------------------------------------
//функция расчёта скорости изменения высоты
T CPUMPSKRotatePropagator::altitude_rate(const T &x,  const T &y,  const T &z,
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



T CPUMPSKRotatePropagator::heightInBaseSystem(const T &x, const T &y, const T &z)
{
  T X = _centerx + x;
  T Y = _centery + y;
  T Z = _centerz + z;
  T dist = RCore::_sqrt(_POW2(X)+_POW2(Y)+_POW2(Z));

  return _heightInMPSK(Y, Z, dist);
}

//приведение величины в интервал [-2*PI, 2*PI]
T bring2interval(const T& val)
{
  return fmod(val, 2.*M_PI);
}




//------------------------------------------------------------------------------
//Прогноз движения объекта в поле Земли на один шаг по методу Моисеева (в МПСК)
CPUBasePropagator::PropagateResult
             CPUMPSKRotatePropagator::propagate(CPUElemParams *iprms,
                                                T &t, const T &dt,
                                                T &x, T &y, T &z,
                                                T &vx, T &vy, T &vz,
                                                T &dvx, T &dvy, T &dvz,
                                                T &pitch, T &raw, T &banking,
                                                T &vpitch, T &vraw, T &vbanking,
                                                T &wx, T &wy, T &wz,
                                                T &dwx, T &dwy, T &dwz,
                                                bool &inatm, T &hght,
                                                const T &atmorange)
{
  if (CPUBasePropagator::propagate(t, dt, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                                   pitch, raw, banking, vpitch, vraw, vbanking,
                                   wx, wy, wz, dwx, dwy, dwz, atmorange)!=
                                   PR_UNIMPLEMENTED)
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }
  if (iprms==nullptr) return PR_INVALID_INPUT_PARAMETERS;

  _params = iprms;


  //метод уточненного расчёта начальной ориентации и скоростей вращения элемента
  if (inatm==false)
  {
     //расчёт тек. высоты
     hght = heightInBaseSystem(x, y, z);
     if (RCore::isNan(hght)) return PR_MATHERROR;
     //скорость изменения высоты
     T ar = altitude_rate(x, y, z, vx, vy, vz);
     if (RCore::isNan(ar)) return PR_MATHERROR;
     //расчёт превышения над граничной высотой
     T hd = hght - _maxheight;

     T dt_ = dt;
     T dtcorr = 0.;
     bool doit = (hd<(-ar*dt_)) && (ar<0.) && (_MFASTABS(hd)>_currency);
     while( doit )
     {
         dtcorr     += -hd/ar;
         //запомнили начальные параметры
         T t1        = t,
           x1        = x,
           y1        = y,
           z1        = z,
           vx1       = vx,
           vy1       = vy,
           vz1       = vz,
           dvx1      = dvx,
           dvy1      = dvy,
           dvz1      = dvz,
           pitch1    = pitch,
           raw1      = raw,
           banking1  = banking,
           vpitch1   = vpitch,
           vraw1     = vraw,
           vbanking1 = vbanking,
           wx1       = wx,
           wy1       = wy,
           wz1       = wz,
           dwx1      = dwx,
           dwy1      = dwy,
           dwz1      = dwz;

        CPUBasePropagator::PropagateResult res =
            __propagate(t1, dtcorr, x1, y1, z1, vx1, vy1, vz1, dvx1, dvy1, dvz1,
                        pitch1, raw1, banking1, vpitch1, vraw1, vbanking1,
                        wx1, wy1, wz1, dwx1, dwy1, dwz1, inatm, hght, atmorange);
        if (res!=PR_OK)
        {
          return res;
        }
        //скорость изменения высоты (altitude rate)
        ar = altitude_rate(x1, y1, z1, vx1, vy1, vz1);
        if (RCore::isNan(ar)) return PR_MATHERROR;
        //расчёт превышения над граничной высотой
        if (RCore::isNan(hght)) hght = heightInBaseSystem(x1, y1, z1);
        hd = hght - _maxheight;

        if ( (ar<0.) && (hd<=0.) )
        {
           //расчёт начальных углов при входе в атмосферу
           _params->calculateInitialAngles(vx1, vy1, vz1, pitch, raw, banking, wx, wy, wz);
           inatm = true;
           bool fl = (hd<(-ar*dt_)) && (ar<0.) && (_MFASTABS(hd)>_currency);
           //процесс уточнения высоты окончен
           if (fl==false)
           {
             dtcorr = dt_ - dtcorr;
            /* t = t1;
             x = x1; y= y1; z = z1; vx = vx1; vy = vy1; vz = vz1;
             dvx = dvx1; dvy = dvy1; dvz = dvz1;
             pitch = pitch1; raw = raw1; banking = banking1;
             vpitch = vpitch1; vraw = vraw1; vbanking = vbanking1;
             wx = wx1; wy = wy1; wz = wz1; dwx = dwx1; dwy = dwy1; dwz = dwz1;*/
             return __propagate(t, dtcorr, x, y, z, vx, vy, vz, dvx, dvy, dvz,
                                pitch, raw, banking, vpitch, vraw, vbanking,
                                wx, wy, wz, dwx, dwy, dwz, inatm, hght, atmorange);
           }
           break;
        }
     }

  }

  CPUBasePropagator::PropagateResult res =
                        __propagate(t, dt, x, y, z, vx, vy, vz, dvx, dvy, dvz,
                                    pitch, raw, banking, vpitch, vraw, vbanking,
                                    wx, wy, wz, dwx, dwy, dwz, inatm, hght, atmorange);
  if (RCore::isNan(hght)) hght = heightInBaseSystem(x, y, z);
  if (res!=PR_OK)
  {
    return res;
  }

  return res;
}

CPUBasePropagator::PropagateResult
           CPUMPSKRotatePropagator::__propagate(T &t, const T &dt,
                                                T &x, T &y, T &z,
                                                T &vx, T &vy, T &vz,
                                                T &dvx, T &dvy, T &dvz,
                                                T &pitch, T &raw, T &banking,
                                                T &vpitch, T &vraw, T &vbanking,
                                                T &wx, T &wy, T &wz,
                                                T &dwx, T &dwy, T &dwz,
                                                bool inatm, T &hght, const T &atmorange)
{
  if (CPUBasePropagator::propagate(t, dt, x,y,z, vx,vy,vz, dvx,dvy,dvz,
                                   pitch, raw, banking, vpitch, vraw, vbanking,
                                   wx, wy, wz, dwx, dwy, dwz, atmorange)!=
                                   PR_UNIMPLEMENTED)
  {
    return PR_INVALID_INPUT_PARAMETERS;
  }
  if (RCore::isNan(hght)) return PR_INVALID_INPUT_PARAMETERS;

  CPUMatrix k0(3, 4),
            k1(3, 4),
            a0(3, 4),
            a1(3, 4),
            lnr(6, 1),
            ang(6, 1);


  T dt2 = dt/2.;


  //расчет линейных и угловых скоростей и ускорений шаг #1
  if (accelerations(x,y,z, vx,vy,vz, dvx,dvy,dvz,
                    pitch, raw, banking, vpitch, vraw, vbanking,
                    wx, wy, wz, dwx, dwy, dwz, inatm, hght, atmorange)==false)
  {
    return PR_MATHERROR;
  }

  k0.v(0,0) = vx;
  k0.v(1,0) = vy;
  k0.v(2,0) = vz;
  k1.v(0,0) = dvx;
  k1.v(1,0) = dvy;
  k1.v(2,0) = dvz;
  lnr[0]    = x + k0.v(0,0)*dt2;
  lnr[1]    = y + k0.v(1,0)*dt2;
  lnr[2]    = z + k0.v(2,0)*dt2;
  lnr[3]    = vx + k1.v(0,0)*dt2;
  lnr[4]    = vy + k1.v(1,0)*dt2;
  lnr[5]    = vz + k1.v(2,0)*dt2;

  a0.v(0,0) = vpitch;
  a0.v(1,0) = vraw;
  a0.v(2,0) = vbanking;
  a1.v(0,0) = dwx;
  a1.v(1,0) = dwy;
  a1.v(2,0) = dwz;
  ang[0]    = pitch   + a0.v(0,0)*dt2;
  ang[1]    = raw     + a0.v(1,0)*dt2;
  ang[2]    = banking + a0.v(2,0)*dt2;
  ang[3]    = wx      + a1.v(0,0)*dt2;
  ang[4]    = wy      + a1.v(1,0)*dt2;
  ang[5]    = wz      + a1.v(2,0)*dt2;


  //шаг #2
  if (accelerations(lnr[0],lnr[1],lnr[2], lnr[3],lnr[4],lnr[5], dvx,dvy,dvz,
                    ang[0], ang[1], ang[2], vpitch, vraw, vbanking,
                    ang[3], ang[4], ang[5], dwx, dwy, dwz, inatm, hght, atmorange)==false)
  {
    return PR_MATHERROR;
  }
  k0.v(0,1) = vx + k1.v(0,0)*dt2;
  k0.v(1,1) = vy + k1.v(1,0)*dt2;
  k0.v(2,1) = vz + k1.v(2,0)*dt2;
  k1.v(0,1) = dvx;
  k1.v(1,1) = dvy;
  k1.v(2,1) = dvz;
  lnr[0]    = x + k0.v(0,1)*dt2;
  lnr[1]    = y + k0.v(1,1)*dt2;
  lnr[2]    = z + k0.v(2,1)*dt2;
  lnr[3]    = vx + k1.v(0,1)*dt2;
  lnr[4]    = vy + k1.v(1,1)*dt2;
  lnr[5]    = vz + k1.v(2,1)*dt2;

  a0.v(0,1) = vpitch;
  a0.v(1,1) = vraw;
  a0.v(2,1) = vbanking;
  a1.v(0,1) = dwx;
  a1.v(1,1) = dwy;
  a1.v(2,1) = dwz;
  ang[0]    = pitch   + a0.v(0,1)*dt2;
  ang[1]    = raw     + a0.v(1,1)*dt2;
  ang[2]    = banking + a0.v(2,1)*dt2;
  ang[3]    = wx      + a1.v(0,1)*dt2;
  ang[4]    = wy      + a1.v(1,1)*dt2;
  ang[5]    = wz      + a1.v(2,1)*dt2;



  //шаг #3
  if (accelerations(lnr[0],lnr[1],lnr[2], lnr[3],lnr[4],lnr[5], dvx,dvy,dvz,
                    ang[0], ang[1], ang[2], vpitch, vraw, vbanking,
                    ang[3], ang[4], ang[5], dwx, dwy, dwz, inatm, hght, atmorange)==false)
  {
    return PR_MATHERROR;
  }
  k0.v(0,2) = vx + k1.v(0,1)*dt2;
  k0.v(1,2) = vy + k1.v(1,1)*dt2;
  k0.v(2,2) = vz + k1.v(2,1)*dt2;
  k1.v(0,2) = dvx;
  k1.v(1,2) = dvy;
  k1.v(2,2) = dvz;
  lnr[0]    = x + k0.v(0,2)*dt;
  lnr[1]    = y + k0.v(1,2)*dt;
  lnr[2]    = z + k0.v(2,2)*dt;
  lnr[3]    = vx + k1.v(0,2)*dt;
  lnr[4]    = vy + k1.v(1,2)*dt;
  lnr[5]    = vz + k1.v(2,2)*dt;

  a0.v(0,2) = vpitch;
  a0.v(1,2) = vraw;
  a0.v(2,2) = vbanking;
  a1.v(0,2) = dwx;
  a1.v(1,2) = dwy;
  a1.v(2,2) = dwz;
  ang[0]    = pitch   + a0.v(0,2)*dt;
  ang[1]    = raw     + a0.v(1,2)*dt;
  ang[2]    = banking + a0.v(2,2)*dt;
  ang[3]    = wx      + a1.v(0,2)*dt;
  ang[4]    = wy      + a1.v(1,2)*dt;
  ang[5]    = wz      + a1.v(2,2)*dt;



  //шаг #4
  if (accelerations(lnr[0],lnr[1],lnr[2], lnr[3],lnr[4],lnr[5], dvx,dvy,dvz,
                    ang[0], ang[1], ang[2], vpitch, vraw, vbanking,
                    ang[3], ang[4], ang[5], dwx, dwy, dwz, inatm, hght, atmorange)==false)
  {
    return PR_MATHERROR;
  }
  k0.v(0,3) = vx + k1.v(0,2)*dt;
  k0.v(1,3) = vy + k1.v(1,2)*dt;
  k0.v(2,3) = vz + k1.v(2,2)*dt;
  k1.v(0,3) = dvx;
  k1.v(1,3) = dvy;
  k1.v(2,3) = dvz;
  t        += dt;
  x        += ( (k0.v(0,0) +2.*k0.v(0,1) + 2.*k0.v(0,2) +  k0.v(0,3))/6.)*dt;
  y        += ( (k0.v(1,0) +2.*k0.v(1,1) + 2.*k0.v(1,2) +  k0.v(1,3))/6.)*dt;
  z        += ( (k0.v(2,0) +2.*k0.v(2,1) + 2.*k0.v(2,2) +  k0.v(2,3))/6.)*dt;
  vx       += ( (k1.v(0,0) +2.*k1.v(0,1) + 2.*k1.v(0,2) +  k1.v(0,3))/6.)*dt;
  vy       += ( (k1.v(1,0) +2.*k1.v(1,1) + 2.*k1.v(1,2) +  k1.v(1,3))/6.)*dt;
  vz       += ( (k1.v(2,0) +2.*k1.v(2,1) + 2.*k1.v(2,2) +  k1.v(2,3))/6.)*dt;

  a0.v(0,3) = vpitch;
  a0.v(1,3) = vraw;
  a0.v(2,3) = vbanking;
  a1.v(0,3) = dwx;
  a1.v(1,3) = dwy;
  a1.v(2,3) = dwz;
  pitch    += ( (a0.v(0,0) +2.*a0.v(0,1) + 2.*a0.v(0,2) +  a0.v(0,3))/6.)*dt;
  raw      += ( (a0.v(1,0) +2.*a0.v(1,1) + 2.*a0.v(1,2) +  a0.v(1,3))/6.)*dt;
  banking  += ( (a0.v(2,0) +2.*a0.v(2,1) + 2.*a0.v(2,2) +  a0.v(2,3))/6.)*dt;
  wx       += ( (a1.v(0,0) +2.*a1.v(0,1) + 2.*a1.v(0,2) +  a1.v(0,3))/6.)*dt;
  wy       += ( (a1.v(1,0) +2.*a1.v(1,1) + 2.*a1.v(1,2) +  a1.v(1,3))/6.)*dt;
  wz       += ( (a1.v(2,0) +2.*a1.v(2,1) + 2.*a1.v(2,2) +  a1.v(2,3))/6.)*dt;

  pitch   = bring2interval(pitch);
  raw     = bring2interval(raw);
  banking = bring2interval(banking);

  if (accelerations(x,y,z, vx,vy,vz, dvx,dvy,dvz,
                    pitch, raw, banking, vpitch, vraw, vbanking,
                    wx, wy, wz, dwx, dwy, dwz, inatm, hght, atmorange)==false)
  {
    return PR_MATHERROR;
  }
  return PR_OK;
}


} //namespace RBallistics

} //namespace RMath

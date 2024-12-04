#include "nan.h"
#include "cpumath.h"
#include "cpumatrix.h"
#include "cpugeodatums.h"

namespace RModels {


namespace RDatums {



//скорость света в вакууме <C> [м\с]
const T lightVelocity = (T)299792458.;
//гравитационная постоянная <f>[м^3\кг*с^2]
const T gravityConst = (T)6.67259e-11;
//средний радиус Земли, [м]
const T mnRadius = (T)6371000.0;

const T consts[6][5] =
{
  //CK42
  {1.08263e-3,  1./298.3,  39861616.79e+7,  6378245.,  0.72921151457e-4},
  //PZ90
  {1.08263e-3,  1./298.257839303,   39860044.e+7,  6378136.,  0.7292115e-4},
  //PZ90_02
  {1.08263e-3, 1./298.257839303, 39860044.e+7, 6378136., 0.7292115e-4},
  //PZ90_11
  {108262.575e-8, 1./298.25784, 39860044.18e+7, 6378136., 0.7292115e-4},
  //GSK2011
  {108262.575e-8, 1./298.2564151, 39860044.15e+7, 6378136.5, 0.7292115e-4},
  //WGS84
  {108263.e-8, 1./298.257223563, 39860050.e+7, 6378137., 0.72921158553e-4}
};

//Датумы
TEarthDatum CK42(TEarthDatum::ET_CK42);
TEarthDatum PZ90(TEarthDatum::ET_PZ90);
TEarthDatum PZ90_02(TEarthDatum::ET_PZ90_02);
TEarthDatum PZ90_11(TEarthDatum::ET_PZ90_11);
TEarthDatum GSK2011(TEarthDatum::ET_GSK2011);
TEarthDatum WGS84(TEarthDatum::ET_WGS84);





TEarthDatum::TEarthDatum(rint32 itp)
{
  if (itp>=ET_USERTYPES) itp = TEarthDatum::ET_CK42;

  _type = itp;
  _J2                       = consts[itp][0];
  _EllipticityEarth         = consts[itp][1];
  _GeocentricalGravityConst = consts[itp][2];
  _BigSemiAxis              = consts[itp][3];
  _AngularVelocityEarth     = consts[itp][4];

  _valid = false;


  //масса Земли(включая массу атмосферы)
  _EarthWeight = _GeocentricalGravityConst/gravityConst;
  T IE = 1./_EllipticityEarth;
  _Nonsphericity = 2.0 / (IE - 1.0) + 1.0 / _POW2(IE - 1.0);

  //производная коэффициента сжатия Земли
  _EllipticityEarthDeriv = (T)2.*_EllipticityEarth -_POW2(_EllipticityEarth);
  //малая полуось, [м]
  _SmallSemiAxis = _BigSemiAxis - _BigSemiAxis*_EllipticityEarth;
  //эксцентриситет
  if (RCore::isEqual(_BigSemiAxis, (T)0.)) return;
  T arg = (T)1. - _POW2(_SmallSemiAxis)/_POW2(_BigSemiAxis);
  if (arg<(T)0.) return;

  _Eccentricity = sqrt(arg);
  //квадрат угловой скорости вращения Земли, [рад^2/с^2]
  _AngularVelocityEarth2 = _POW2(_AngularVelocityEarth);

  //Параметры для учёта влияния гравитационного поля Земли
  const T cz = 0.5*_J2*_POW2(_BigSemiAxis);
  _a00   = _GeocentricalGravityConst/mnRadius;
  _a20   = -2.*cz*_GeocentricalGravityConst/_POW3(mnRadius);

  _valid = true;
}
TEarthDatum::TEarthDatum(rint32 itp, const T &ij2, const T &iea,
                         const T &iggc, const T &ibsa, const T &iave)
                         :_GeocentricalGravityConst(iggc),_EllipticityEarth(iea),
                       _BigSemiAxis(ibsa),_AngularVelocityEarth(iave),_J2(ij2)
{
  _valid = false;

  _type = itp;

  if (itp<ET_USERTYPES)
  {
    _J2                       = consts[itp][0];
    _EllipticityEarth         = consts[itp][1];
    _GeocentricalGravityConst = consts[itp][2];
    _BigSemiAxis              = consts[itp][3];
    _AngularVelocityEarth     = consts[itp][4];
  }

  //масса Земли(включая массу атмосферы)
  _EarthWeight = _GeocentricalGravityConst/gravityConst;
  //производная коэффициента сжатия Земли
  _EllipticityEarthDeriv = (T)2.*_EllipticityEarth -_POW2(_EllipticityEarth);
  //малая полуось, [м]
  _SmallSemiAxis = _BigSemiAxis - _BigSemiAxis*_EllipticityEarth;
  //эксцентриситет
  if (RCore::isEqual(_BigSemiAxis, (T)0.)) return;
  T arg = (T)1. - _POW2(_SmallSemiAxis)/_POW2(_BigSemiAxis);
  if (arg<(T)0.) return;

  _Eccentricity = sqrt(arg);
  T IE = 1./iea;
  _Nonsphericity = 2.0 / (IE - 1.0) + 1.0 / _POW2(IE - 1.0);
  //квадрат угловой скорости вращения Земли, [рад^2/с^2]
  _AngularVelocityEarth2 = _POW2(_AngularVelocityEarth);

  //Параметры для учёта влияния гравитационного поля Земли
  const T cz = 0.5*_J2*_POW2(_BigSemiAxis);
  _a00   = _GeocentricalGravityConst/mnRadius;
  _a20   = -2.*cz*_GeocentricalGravityConst/_POW3(mnRadius);

  _valid = true;
}

TEarthDatum::TEarthDatum(const TEarthDatum& copy)
{
  _valid = false;

  _GeocentricalGravityConst = 0.;
  _EllipticityEarth         = 0.;
  _BigSemiAxis              = 0.;
  _AngularVelocityEarth     = 0.;
  _J2                       = 0.;
  _type                     = 0;
  _EarthWeight              = 0.;
  _EllipticityEarthDeriv    = 0.;
  _SmallSemiAxis            = 0.;
  _Eccentricity             = 0.;
  _Nonsphericity            = 0.;
  _a00                      = 0.;
  _a20                      = 0.;
  _AngularVelocityEarth2    = 0.;

  *this      = copy;
}

//------------------------------------------------------------------------------
//копирование датума (перегрузка оператора присваивания)                      +
TEarthDatum & TEarthDatum::operator=(const TEarthDatum &copy)
{
  if (this == &copy) return *this;

  _valid                    = copy._valid;
  _GeocentricalGravityConst = copy._GeocentricalGravityConst;
  _a00                      = copy._a00;
  _a20                      = copy._a20;
  _EllipticityEarth         = copy._EllipticityEarth;
  _BigSemiAxis              = copy._BigSemiAxis;
  _AngularVelocityEarth     = copy._AngularVelocityEarth;
  _J2                       = copy._J2;
  _type                     = copy._type;
  _EarthWeight              = copy._EarthWeight;
  _EllipticityEarthDeriv    = copy._EllipticityEarthDeriv;
  _SmallSemiAxis            = copy._SmallSemiAxis;
  _Eccentricity             = copy._Eccentricity;
  _Nonsphericity            = copy._Nonsphericity;
  _AngularVelocityEarth2    = copy._AngularVelocityEarth2;

  return *this;
}


//------------------------------------------------------------------------------
//корректно ли создан объект
bool TEarthDatum::isValid()const
{
    return _valid;
}


//------------------------------------------------------------------------------
//тип
rint32 TEarthDatum::type()const
{
  return _type;
}
//------------------------------------------------------------------------------
//коэффициент сжатия Земли
T TEarthDatum::EllipticityEarth()const
{
  return _EllipticityEarth;
}
//------------------------------------------------------------------------------
//производная коэффициента сжатия Земли
T TEarthDatum::EllipticityEarthDeriv()const
{
  return _EllipticityEarthDeriv;
}
//------------------------------------------------------------------------------
//большая полуось(экваториальный радиус)
T TEarthDatum::BigSemiAxis()const
{
  return _BigSemiAxis;
}
//------------------------------------------------------------------------------
//малая полуось
T TEarthDatum::SmallSemiAxis()const
{
  return _SmallSemiAxis;
}
//------------------------------------------------------------------------------
//эксцентриситет
T TEarthDatum::Eccentricity()const
{
  return _Eccentricity;
}
//------------------------------------------------------------------------------
//Коэффициент несферичности Земли
T TEarthDatum::Nonsphericity()const
{
  return _Nonsphericity;
}
//------------------------------------------------------------------------------
//геоцентрическая гравитационная постоянная Земли
T TEarthDatum::GravityConstant()const
{
    return _GeocentricalGravityConst;
}

//------------------------------------------------------------------------------
//Параметры для учёта влияния гравитационного поля Земли
T TEarthDatum::a00() const
{
  return _a00;
}
T TEarthDatum::a20() const
{
  return _a20;
}

//------------------------------------------------------------------------------
//угловая скорость вращения Земли
T TEarthDatum::AngularVelocity()const
{
  return _AngularVelocityEarth;
}
//------------------------------------------------------------------------------
//квадрат угловой скорости вращения Земли
T TEarthDatum::AngularVelocity2()const
{
  return _AngularVelocityEarth2;
}

//------------------------------------------------------------------------------
//средний радиус Земли
T TEarthDatum::meanRadius()
{
  return mnRadius;
}

//------------------------------------------------------------------------------
//масса Земли(включая массу атмосферы)
T TEarthDatum::EarthWeight()const
{
  return _EarthWeight;
}
//------------------------------------------------------------------------------
//коэффициент 2-й зональной гармоники
T TEarthDatum::J2()const
{
  return _J2;
}

//------------------------------------------------------------------------------
//отношение массы Земли к массе Луны.
T TEarthDatum::MoonMassKoef()
{
  return 82.300587;
}
//------------------------------------------------------------------------------
//ускорение свободного падения на Земле
T TEarthDatum::G0()
{
  return 9.80665;
}







//трансформанты координат между ПЗ
T transformants[10][7] =
{
   // m,  dX, dY, dZ, Wx, Wy, Wz
   //+СК-42 ==> ПЗ-90
   {  0.,     25.,    -141.,  -80.,  0.,   _MDEG2RAD(-0.35/3600.),   _MDEG2RAD(-0.66/3600.)},
   //+СК-42 ==> ПЗ-90.02
   {  -0.22e-6,   23.93,     -141.03,  -79.98,  0.,   _MDEG2RAD(-0.35/3600.),   _MDEG2RAD(-0.79/3600.)},
   //+СК-42 ==> ПЗ-90.11
   {  -0.228e-6,     23.557,    -140.844,       -79.778,  _MDEG2RAD(-0.0023/3600.),   _MDEG2RAD(-0.34646/3600.),  _MDEG2RAD(-0.79421/3600.)},

   //+ПЗ-90 ==> ПЗ-90.02
   {  -0.22e-6,   -1.07,     -0.03,      0.02,    0.,  0.,   _MDEG2RAD(-0.13/3600.) },
   //+ПЗ-90 ==> ПЗ-90.11
   {  -0.228e-6,   -1.443,    0.156,     0.222,   _MDEG2RAD(-0.0023/3600.),  _MDEG2RAD(0.00354/3600.),   _MDEG2RAD(-0.13421/3600.) },

   //+ПЗ-90.02 ==> ПЗ-90.11
   {  -0.008e-6, -0.373,     0.186,    0.202,    _MDEG2RAD(-0.0023/3600.),   _MDEG2RAD(0.00354/3600.),   _MDEG2RAD(-0.00421/3600.) },

   //+ГСК-2011 ==> ПЗ-90.11
   { -0.0006e-6,   0.,   0.014,  -0.008,    _MDEG2RAD(-0.000562/3600.),  _MDEG2RAD(-0.000019/3600.),   _MDEG2RAD(0.000053/3600.)  },

   //?WGS-84 ==> ПЗ-90
   {  -0.12e-6,  -1.1,  -0.3,  -0.9,      0.,  0.,  _MDEG2RAD(-0.02/3600.) },
   //+WGS-84 ==> ПЗ-90.02
   {  0.,  0.36,  -0.08,  -0.18,      0.,  0.,  0. },
   //+WGS-84 ==> ПЗ-90.11
   {  -0.008e-6,  -0.013,  0.106,   0.022,  _MDEG2RAD(-0.0023/3600.),   _MDEG2RAD(0.00354/3600.),  _MDEG2RAD(-0.00421/3600.) }

};

//------------------------------------------------------------------------------
//Трансформация координат между параметрами Земли
bool TEarthDatum::transformByEarthDatums(T &ix, T &iy, T &iz,
                                         const TEarthDatum::CEarthTypes ifrom,
                                         const TEarthDatum::CEarthTypes ito)
{
  if (RCore::isNan(ix) || RCore::isNan(iy) || RCore::isNan(iz)) return false;
  if (ifrom==ito) return true; //преобразование не требуется

  T sign = 1.;
  rint32 idx = -1;

  if (idx<0) return false;

  switch (ifrom)
  {
    case ET_CK42:
    {
      if (ito==ET_PZ90)
      {
        idx = 0;
      }
      else if (ito==ET_PZ90_02)
      {
        idx = 1;
      }
      else if (ito==ET_PZ90_11)
      {
        idx = 2;
      }
//      else if (ito==ET_WGS84)
//      {
//        idx = -1;
//      }
//      else if (ito==ET_GSK2011)
//      {
//        idx = -1;
//      }
      break;
    }
    case ET_PZ90:
    {
      if (ito==ET_CK42)
      {
        idx  = 0;
        sign = -1.;
      }
      else if (ito==ET_PZ90_02)
      {
        idx = 3;
      }
      else if (ito==ET_PZ90_11)
      {
        idx = 4;
      }
      else if (ito==ET_WGS84)
      {
        idx = 7;
        sign = -1.;
      }
//      else if (ito==ET_GSK2011)
//      {
//        idx = -1;
//      }
      break;
    }
    case ET_PZ90_02:
    {
      if (ito==ET_CK42)
      {
        idx  = 1;
        sign = -1.;
      }
      else if (ito==ET_PZ90)
      {
        sign = -1.;
        idx = 3;
      }
      else if (ito==ET_PZ90_11)
      {
        idx = 5;
      }
      else if (ito==ET_WGS84)
      {
        idx = 8;
        sign = -1.;
      }
//      else if (ito==ET_GSK2011)
//      {
//        idx = -1;
//      }
      break;
    }
    case ET_PZ90_11:
    {
      if (ito==ET_CK42)
      {
        idx  = 2;
        sign = -1.;
      }
      else if (ito==ET_PZ90)
      {
        sign = -1.;
        idx = 4;
      }
      else if (ito==ET_PZ90_02)
      {
        sign = -1.;
        idx = 5;
      }
      else if (ito==ET_WGS84)
      {
        sign = -1.;
        idx = 9;
      }
      else if (ito==ET_GSK2011)
      {
        sign = -1.;
        idx  = 6;
      }
      break;
    }
    case ET_GSK2011:
    {
      if (ito==ET_PZ90_11)
      {
        idx  = 6;
      }
//      else if (ito==ET_PZ90)
//      {
//      }
//      else if (ito==ET_PZ90_02)
//      {
//      }
//      else if (ito==ET_CK42)
//      {
//      }
//      else if (ito==ET_WGS84)
//      {
//      }
      break;
    }
    case ET_WGS84:
    {
      if (ito==ET_PZ90)
      {
        idx = 7;
      }
      else if (ito==ET_PZ90_02)
      {
        idx = 8;
      }
      else if (ito==ET_PZ90_11)
      {
        idx = 9;
      }
//      else if (ito==ET_CK42)
//      {
//      }
//      else if (ito==ET_GSK2011)
//      {
//      }
      break;
    }

    default:idx = -1;
  }
  //масштабирование
  T mult = (1. + sign*transformants[idx][0]);
  //смещение
  RMatrix::CPUMatrix shift(1, 3);
  shift[0] = sign*transformants[idx][1];
  shift[1] = sign*transformants[idx][2];
  shift[2] = sign*transformants[idx][3];
  //поворот
  RMatrix::CPUMatrix rotor(3, 3, true);
  rotor.v(2,1) = sign*transformants[idx][4];
  rotor.v(1,2) = -sign*transformants[idx][4];
  rotor.v(2,0) = -sign*transformants[idx][5];
  rotor.v(0,2) = sign*transformants[idx][5];
  rotor.v(1,0) = sign*transformants[idx][6];
  rotor.v(0,1) = -sign*transformants[idx][6];
  //ИД
  RMatrix::CPUMatrix pos(1, 3);
  pos[0] = ix;
  pos[1] = iy;
  pos[2] = iz;

  pos = mult*rotor*pos + shift;
  ix = pos[0];
  iy = pos[1];
  iz = pos[2];

  return true;
}






} //namespace RDatums

} //namespace RModels

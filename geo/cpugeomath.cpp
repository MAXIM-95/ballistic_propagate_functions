#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cpumath.h"
#include "cpumatrix.h"
#include "cpugeomath.h"
#include "cpugeodatums.h"

namespace RMath {



namespace RGeoMath {




//------------------------------------------------------------------------------
//                ФУНКЦИИ ПРЕОБРАЗОВАНИЯ СИСТЕМ КООРДИНАТ                      |
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Преобразование Геодезических координат в Геоцентрические                     |
//вх.данные  [b,l,h] геодезические (широта,[рад]; долгота,[рад]; высота,[м])   |
//            earth - параметры земли (по умолчанию СК-42)                     |
//вых.данные [Bgc,Lgc,Hgc] геодезические (широта,[рад]; долгота,[рад]; выс.,[м]|
bool geodez2geocentrical(const T &B, const T &L, const T &H,
                         T &Bgc, T &Lgc, T &Hgc, const TEarthDatum &earth)
{
#ifdef G_ASSERT_SHEME
  assert(earth.SmallSemiAxis()!=(T)0.);
#else
  if (earth.SmallSemiAxis()==(T)0.) return false;
#endif

  Bgc = RCore::_atan(_POW2(earth.BigSemiAxis())/_POW2(earth.SmallSemiAxis())*
        RCore::_tan(B));
  Lgc = L;
  Hgc = H;

  return true;
}

bool geodez2geocentrical(const T &B, const T &L, const T &H,
                         const T &dB, const T &dL, const T &dH,
                         T &Bgc, T &Lgc, T &Hgc,
                         T &dBgc, T &dLgc, T &dHgc, const TEarthDatum &earth)
{
  if (!geodez2geocentrical(B, L, H, Bgc, Lgc, Hgc, earth)) return false;

  //расчёт производных
  if (RCore::isNan(dB) || RCore::isNan(dL) || RCore::isNan(dH))
  {
    dBgc = dLgc = dHgc = NAN;
    return false;
  }

  bool res = geodez2geocentrical(B+dB, dL, dH, dBgc, dLgc, dHgc,  earth);
  if (res)
  {
     dBgc -= Bgc;
  }
  return res;
}


//------------------------------------------------------------------------------
//Преобразование Геоцентрических координат в Геодезические                     |
//вх.данные  [Bgc,Lgc,Hgc] геоцентрические (широта,[рад]; долгота,[рад];выс,[м]|
//[dBgc,dLgc,dHgc] геоц. произв.(широты,[рад/c]; долготы,[рад/c]; высоты,[м/c])|
//            earth - параметры земли (по умолчанию СК-42)                     |
//вых.данные [B,L,H] геодезические (широта,[рад]; долгота,[рад]; высота,[м])   |
//[dB, dL, dH]геодез. производные(широты,[рад/c]; долготы,[рад/c]; высоты,[м/c])
bool geocentrical2geodez(const T &Bgc, const T &Lgc, const T &Hgc,
                         T &B, T &L, T &H, const TEarthDatum &earth)
{
#ifdef G_ASSERT_SHEME
  assert(earth.BigSemiAxis()!=(T)0.);
#else
  if (earth.BigSemiAxis()==(T)0.) return false;
#endif

  B = atan(_POW2(earth.SmallSemiAxis())/_POW2(earth.BigSemiAxis())*tan(Bgc));
  L = Lgc;
  H = Hgc;

  return true;
}

bool geocentrical2geodez(const T &Bgc, const T &Lgc, const T &Hgc,
                         const T &dBgc, const T &dLgc, const T &dHgc,
                         T &B, T &L, T &H, T &dB, T &dL, T &dH,
                         const TEarthDatum &earth)
{
  if (!geocentrical2geodez(Bgc, Lgc, Hgc, B, L, H, earth)) return false;

  //расчёт производных
  if (RCore::isNan(dBgc) || RCore::isNan(dLgc) || RCore::isNan(dHgc))
  {
    dB = dL = dH = NAN;
    return false;
  }

  bool res = geocentrical2geodez(dBgc + Bgc, dLgc, dHgc, dB, dL, dH, earth);
  if (res)
  {
     dB -= B;
  }
  return res;
}



//------------------------------------------------------------------------------
//Преобразование из Геодезической СК в ПГСК                                    |
//вх.данные  [b,l,h] (широта,[рад]; долгота,[рад]; высота,[м])                 |
//        [db,dl,dh] скорости по (широте,[рад\с]; долгота,[рад\с]; высоте,[м\с])
//            earth - параметры земли (по умолчанию СК-42)                     |
//вых.данные [x,y,z] (координаты в ПГСК, [м])                                  |
bool geodez2pgsk(const T &B, const T &L, const T &H,
                 T &X, T &Y, T &Z,  const TEarthDatum &earth)
{
  T  SINL = RCore::_sin(L), COSL = RCore::_cos(L),
     SINB = RCore::_sin(B), COSB = RCore::_cos(B),
     e2   = earth.EllipticityEarthDeriv(),
     tmp  = RCore::_sqrt((T)1.-e2*_POW2(SINB));

#ifdef G_ASSERT_SHEME
  assert(tmp!=(T)0.);
#else
  if (tmp==(T)0.) return false;
#endif

  T  M    = earth.BigSemiAxis()/tmp;

  X = (M + H)*COSB*COSL;
  Y = (M + H)*COSB*SINL;
  Z = (M*((T)1. - e2) + H)*SINB;

  return true;
}

bool geodez2pgsk(const T &B, const T &L, const T &H,
                 const T &dB, const T &dL, const T &dH,
                 T &X, T &Y, T &Z, T &VX, T &VY, T &VZ,
                 const TEarthDatum &earth)
{
  if (!geodez2pgsk(B, L, H, X, Y, Z, earth)) return false;

  //расчёт производных
  if (RCore::isNan(dB) || RCore::isNan(dL) || RCore::isNan(dH))
  {
    VX = VY = VZ = NAN;
    return false;
  }

  bool res = geodez2pgsk(B+dB, L+dL, H+dH, VX, VY, VZ,  earth);
  if (res)
  {
     VX -= X;
     VY -= Y;
     VZ -= Z;
  }
  return res;
}



//------------------------------------------------------------------------------
//Преобразование из ПГСК в Геодезическую СК                                    |
//вх.данные [x,y,z] (координаты в ПГСК, [м])                                   |
//           earth - параметры земли (по умолчанию СК-42)                      |
//вых.данные  [b,l,h] (широта,[рад]; долгота,[рад]; высота,[м])                |
bool pgsk2geodez(const T &X, const T &Y, const T &Z,
                 T &B, T &L, T &H, const TEarthDatum &earth)
{
   T r2 = _POW2(X) + _POW2(Y),
     b  = earth.SmallSemiAxis(),
     a  = earth.BigSemiAxis();
#ifdef G_ASSERT_SHEME
  assert(a!=0.0);
#else
  if (a==(T)0.) return false;
#endif
   if (r2==0.0)
   {
      B = (Z>=0.0)? M_PI_2 : -M_PI_2;
      L = 0.0;
      H = _MFASTABS(Z) - b;

      return true;
   }

   const T r = RCore::_sqrt(r2);
   L = RCore::_atan2(Y, X);
   if (Z < 0.0) b = -b;
   r2 = b/a;
#ifdef G_ASSERT_SHEME
   assert(r!=0.0);
#else
   if (r==(T)0.) return false;
#endif
   const T TT  = ((Z + b)*r2 - a)/r,
           TT2 = _POW2(TT),
           F   = ((Z - b)*r2 + a)/r,
           P   = (TT*F + 1.0)*4.0/3.0,
           Q   = (TT2 - _POW2(F))*2.0,
           D   = _POW3(P) + _POW2(Q);

    T  v;
    if (D >= 0.0)
    {
        T s = RCore::_sqrt(D) + Q;
        v = RCore::_exp(RCore::_log(_MFASTABS(s))*(1.0/3.0));
        s = (s < 0.0) ?-v :v;
#ifdef G_ASSERT_SHEME
   assert(s!=0.0);
#else
   if (s==(T)0.) return false;
#endif
        v = P/s - s;
        v = -(Q + Q + _POW3(v))/(3.0*P);
    }
    else
    {
#ifdef G_ASSERT_SHEME
        assert((-P)>=0.0);
#else
        if ((-P)<0.0) return false;
#endif
       const T sqrtp = RCore::_sqrt(-P);
       v = (P*sqrtp);
#ifdef G_ASSERT_SHEME
       assert(v!=0.0);
#else
       if (v==(T)0.) return false;
#endif
       v = Q/v *3.0;
#ifdef G_ASSERT_SHEME
       assert( (v >= -1.0) && (v <= 1.0) );
#else
       if ( (v < -1.0) || (v > 1.0) ) return false;
#endif
       v = 2.0*sqrtp*RCore::_cos(RCore::_acos(v));
    }
    T   G = 0.5*(TT + RCore::_sqrt(TT2 + v)),
        t = RCore::_sqrt(_POW2(G) + (F - v*G)/(G + G - TT)) - G;

    B = RCore::_atan((1.0 - _POW2(t))*a/(2.0*b*t));
    T sn = RCore::_sin(B),
      cs = RCore::_cos(B);
    H = (r - a*t)*cs + (Z - b)*sn;

    return true;
}


bool pgsk2geodez(const T &X, const T &Y, const T &Z,
                     const T &VX, const T &VY, const T &VZ,
                     T &B, T &L, T &H, T &dB, T &dL, T &dH,
                     const TEarthDatum &earth)
{
  if (!pgsk2geodez(X, Y, Z, B, L, H, earth)) return false;

  //расчёт производных
  if (RCore::isNan(VX) || RCore::isNan(VY) || RCore::isNan(VZ))
  {
    dB = dL = dH = NAN;
    return false;
  }

  bool res = pgsk2geodez(X+VX, Y+VY, Z+VZ, dB, dL, dH, earth);
  if (res)
  {
     dB -= B;
     dL -= L;
     dH -= H;
  }
  return res;
}




//------------------------------------------------------------------------------
//Преобразование из РТСК в МПСК                                                |
//вх.данные [R,B,E,dR,dB,dE] (дальность, [м]; азимут, [рад]; угол места, [рад])|
//(сферические РТСК координаты и их производные)                               |
//вых.данные  [X,Y,Z,VX,VY,VZ] прямоугольные МПСК координаты и их производные  |
bool rtsk2mpsk(const T &R,  const T &B,  const T &E,
               const T &dR, const T &dB, const T &dE,
               T &X, T &Y, T &Z, T &VX, T &VY, T &VZ)
{
  T  SINE = RCore::_sin(E), COSE = RCore::_cos(E),
     SINB = RCore::_sin(B), COSB = RCore::_cos(B);

  X  = R*COSE*COSB;
  Y  = R*SINE;
  Z  = R*SINB*COSE;

  //расчёт производных
  if ( (RCore::isNan(dR)==false) && (RCore::isNan(dB)==false) &&
       (RCore::isNan(dE)==false) )
  {
    VX = dR*COSB*COSE - dB*Z - Y*COSB*dE;
    VY = dR*SINE + dE*R*COSE;
    VZ = dR*SINB*COSE + dB*X - Y*SINB*dE;
  }
  else
  {
    VX = VY = VZ = NAN;
  }

  return true;
}



//------------------------------------------------------------------------------
//Преобразование из МПСК в РТСК                                                |
//вх.данные [X,Y,Z,VX,VY,VZ] прямоуг. МПСК координаты [м] и их производные[м/с]|
//вых.данные [R,B,E,dR,dB,dE] (дальность,[м]; азимут, [рад]; угол места, [рад])|
//(сферические РТСК координаты и их производные)                               |
bool mpsk2rtsk(const T &X, const T &Y, const T &Z,
               const T &VX, const T &VY, const T &VZ,
               T &R,  T &B,  T &E, T &dR, T &dB, T &dE)
{
  T aa = X*X + Z*Z;
  R = RCore::_sqrt(aa + Y*Y);
  if (R<_MFASTABS(Y)) R = _MFASTABS(Y);

  bool nanfl = RCore::isNan(VX) || RCore::isNan(VY) || RCore::isNan(VZ);
  if (nanfl)
  {
    dR = dB = dE = NAN;
  }

  if (RCore::isEqual(R, (T)0.))
  {
    E = 0.;
    if (!nanfl) dR = 0.;
  }
  else
  {
    E  = RCore::_asin(Y/R);
    if (!nanfl) dR = (X*VX + Y*VY + Z*VZ)/R;
  }

  if (fabs(aa)>0.1e-12)
  {
    B  = RCore::atan12(Z, X);
    if (!nanfl) dB = (X*VZ - Z*VX)/aa;
    if (!nanfl) dE = (RCore::isEqual(B, (T)0.)) ?0.
                                                :(VY - Y*dR/R)/RCore::_sqrt(aa);
  }
  else
  {
    B  = RCore::atan12(VZ, VX);
    if (!nanfl) dB = 0.;
    if (!nanfl) dE = (RCore::isEqual(R, (T)0.))  ?0.
                                          :(RCore::_cos(E)*VY - RCore::_sin(E) *
                                     (RCore::_cos(B)*VX + RCore::_sin(B)*VZ))/R;
  }

  return true;
}





//------------------------------------------------------------------------------
//расчёт матрицы поворота для перехода от ПГСК к МПСК                          |
//вх.данные [B,L,Q] геодезические широта, долгота [рад], и азимут РСН [рад]    |
//вых.данные[e11,21,e31,e12,e22,e32,e13,e23,e33] элементы матрицы поворота     |
bool pgsk2mpskRotationMatrix(const T &B, const T &L, const T &Q,
                       T &e11, T &e21, T &e31,
                       T &e12, T &e22, T &e32,
                       T &e13, T &e23, T &e33)
{
  T  SINL = RCore::_sin(L), COSL = RCore::_cos(L),
     SINB = RCore::_sin(B), COSB = RCore::_cos(B),
     SINQ = RCore::_sin(Q), COSQ = RCore::_cos(Q);

  e11 = -SINL*SINQ-COSL*SINB*COSQ;
  e21 = COSL*COSB;
  e31 = -SINL*COSQ+COSL*SINB*SINQ;
  e12 = COSL*SINQ-SINL*SINB*COSQ;
  e22 = SINL*COSB;
  e32 = COSL*COSQ+SINL*SINB*SINQ;
  e13 = COSB*COSQ;
  e23 = SINB;
  e33 = -COSB*SINQ;

  return true;
}


//------------------------------------------------------------------------------
//расчёт матрицы перехода от ПГСК к МПСК                                       |
//вх.данные [B,L,H, Q] геодезические широта, долгота [рад],                    |
//                     и высота центра МПСК, азимут РСН [рад]                  |
//вых.данные[e11,21,e31,e12,e22,e32,e13,e23,e33] элементы матрицы перехода     |
//          [shX, shY, shZ] смещение центра координат МПСК                     |
bool pgsk2mpskMatr(const T &B, const T &L, const T &H, const T &Q,
                   T &e11, T &e21, T &e31, T &e12, T &e22, T &e32,
                   T &e13, T &e23, T &e33, T &shX, T &shY, T &shZ,
                   const TEarthDatum &earth)
{
  if (!geodez2pgsk(B, L, H, shX, shY, shZ, earth)) return false;

  return pgsk2mpskRotationMatrix(B, L, Q,
                                 e11, e21, e31, e12, e22, e32, e13, e23, e33);
}




//------------------------------------------------------------------------------
//Преобразование из ПГСК в МПСК                                                |
//вх.данные: [shX, shY, shZ] смещение центра координат МПСК                    |
//           [e11,21,e31,e12,e22,e32,e13,e23,e33] элементы матрицы перехода    |
//           [iX,iY,iZ,iVX,iVY,iVZ]  ПГСК координаты [м] и их производные[м/с] |
//вых.данные: [oX,oY,oZ,oVX,oVY,oVZ] МПСК координаты [м] и их производные[м/с] |
bool pgsk2mpsk(const T &iX, const T &iY, const T &iZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ)
{
   //сдвиг СК в точку стояния МПСК
   T tX = iX - shX,
     tY = iY - shY,
     tZ = iZ - shZ;

   //расчёт координат положения
   oX = tX*e11 + tY*e12 + tZ*e13;
   oY = tX*e21 + tY*e22 + tZ*e23;
   oZ = tX*e31 + tY*e32 + tZ*e33;

   return true;
}
bool pgsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ)
{
   //сдвиг СК в точку стояния МПСК
   T tX = iX - shX,
     tY = iY - shY,
     tZ = iZ - shZ;

   //расчёт координат положения
   oX = tX*e11 + tY*e12 + tZ*e13;
   oY = tX*e21 + tY*e22 + tZ*e23;
   oZ = tX*e31 + tY*e32 + tZ*e33;

   //расчёт скоростей
   if ( (RCore::isNan(iVX)==false) && (RCore::isNan(iVY)==false) &&
        (RCore::isNan(iVZ)==false) )
   {
     oVX = iVX*e11 + iVY*e12 + iVZ*e13;
     oVY = iVX*e21 + iVY*e22 + iVZ*e23;
     oVZ = iVX*e31 + iVY*e32 + iVZ*e33;
   }
   else
   {
     oVX = oVY = oVZ = NAN;
   }

   return true;
}


//------------------------------------------------------------------------------
//Преобразование из МПСК в ПГСК                                                |
//вх.данные: [shX, shY, shZ] смещение центра координат МПСК                    |
//           [e11,21,e31,e12,e22,e32,e13,e23,e33] элементы матрицы перехода    |
//           [iX,iY,iZ,iVX,iVY,iVZ]  ПГСК координаты [м] и их производные[м/с] |
//вых.данные: [oX,oY,oZ,oVX,oVY,oVZ] МПСК координаты [м] и их производные[м/с] |
bool mpsk2pgsk(const T &iX, const T &iY, const T &iZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ)
{
   //расчёт координат положения
   oX = iX*e11 + iY*e21 + iZ*e31;
   oY = iX*e12 + iY*e22 + iZ*e32;
   oZ = iX*e13 + iY*e23 + iZ*e33;

   //сдвиг СК в точку стояния МПСК
   oX += shX; oY += shY; oZ += shZ;

   return true;
}
bool mpsk2pgsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ)
{
   //расчёт координат положения
   oX = iX*e11 + iY*e21 + iZ*e31;
   oY = iX*e12 + iY*e22 + iZ*e32;
   oZ = iX*e13 + iY*e23 + iZ*e33;

   //сдвиг СК в точку стояния МПСК
   oX += shX; oY += shY; oZ += shZ;

   //расчёт скоростей
   if ( (RCore::isNan(iVX)==false) && (RCore::isNan(iVY)==false) &&
        (RCore::isNan(iVZ)==false) )
   {
     oVX = iVX*e11 + iVY*e21 + iVZ*e31;
     oVY = iVX*e12 + iVY*e22 + iVZ*e32;
     oVZ = iVX*e13 + iVY*e23 + iVZ*e33;
   }
   else
   {
     oVX = oVY = oVZ = NAN;
   }

   return true;
}






//------------------------------------------------------------------------------
//Преобразование координат из РТСК в АПСК                                      |
// вх.данные:     [r[м],b[рад],e[рад],dR[м\с],dB[рад\с],dE[рад\с]]             |
// вых.данные:    [x[м],y[м],z[м],VX[м\с],VY[м\с],VZ[м\с]]                     |
bool rtsk2apsk(const T &r, const T &b, const T &e,
               const T &dR, const T &dB, const T &dE,
               T &x, T &y, T &z, T &VX, T &VY, T &VZ)
{
   const T c1 = RCore::_cos(b), c2 = RCore::_sin(b),
           c3 = RCore::_cos(e), c4 = RCore::_sin(e),
           c5 = -c2*c3, c6 = c1*c3, c7 = c2*c4, c8 = c1*c4,
           c9 = -dB*c6 + dE*c7, c10 = dB*c5 - dE*c8, c11 = dE*c3;

   x  = r*c5;
   y  = r*c4;
   z  = r*c6;

   //расчёт производных
   if ( (RCore::isNan(dR)==false) && (RCore::isNan(dB)==false) &&
        (RCore::isNan(dE)==false) )
   {
     VX = dR*c5 + r*c9;
     VY = dR*c4 + r*c11;
     VZ = dR*c6 + r*c10;
   }
   else
   {
     VX = VY = VZ = NAN;
   }

   return true;
}


//------------------------------------------------------------------------------
//Преобразование координат из АПСК в РТСК                                      |
// вх.данные:   [ x[м],y[м],z[м],VX[м\с],VY[м\с],VZ[м\с] ]                     |
// вых.данные:  [ r[м], b[рад], e[рад], dR[м\с], dB[рад\с], dE[рад\с] ]        |
bool apsk2rtsk(const T &x,  const T &y,  const T &z,
               const T &VX, const T &VY, const T &VZ,
               T &r, T &b, T &e, T &dR, T &dB, T &dE)
{
   r = RCore::_sqrt(x*x + y*y + z*z);
   b = 0.;
   e = 0.;

   bool nanfl = RCore::isNan(VX) || RCore::isNan(VY) || RCore::isNan(VZ);
   if (nanfl)
   {
     dR = dB = dE = NAN;
   }
   else
   {
     dR = dB = dE = 0.;
   }

   if (RCore::isEqual(r, (T)0.)) return false;

   T arg = y/r;

   if (arg<=(T)1. && arg>=(T)-1.) e = RCore::_asin(arg);

   const T tmp = r * RCore::_cos(e);
   if (RCore::isEqual(tmp, (T)0.)) return false;

   arg = z/tmp;

   if (arg<=(T)1. && arg>=(T)-1.)
   {
     if (x<=(T)0.) b = RCore::_acos(arg);
     else b = (T)2.*(T)M_PI - RCore::_acos(arg);
   }

   if (!nanfl)
   {
     dR = (x*VX + y*VY + z*VZ)/r;
     dB = (-VX*RCore::_cos(b) - VZ*RCore::_sin(b))/tmp;
     dE = (VY - dR*RCore::_sin(e))/tmp;
   }

   return true;
}




//------------------------------------------------------------------------------
//Определение дальности по Земле от точки X до точки Y                         |
//вх.данные: [Bx, Lx, Hx] - геодезические координаты точки Х (радианы и метры) |
//           [By, Ly, Hy] - геодезические координаты точки Y                   |
//вых.данные: дальности по Земле от точки X до точки Y, [м]                    |
T sphereDistance(const T &Bx, const T &Lx, const T &Hx,
                 const T &By, const T &Ly, const T &Hy, const TEarthDatum &earth)
{
   if (RCore::isEqual(earth.BigSemiAxis(), (T)0.)) return false;

   T Bx1, Lx1, Hx1;
   geocentrical2geodez(Bx, Lx, Hx, Bx1, Lx1, Hx1, earth);
   T By1, Ly1, Hy1;
   geocentrical2geodez(By, Ly, Hy, By1, Ly1, Hy1, earth);

   T  x1, y1, z1, x2, y2, z2;

   if (!geodez2pgsk(Bx1, Lx1, Hx1, x1, y1, z1, earth) ||
       !geodez2pgsk(By1, Ly1, Hy1, x2, y2, z2, earth)) return NAN;

   //нормирование координат делением на большую полуось земного эллипсоида
   T norm = (T)1.0/earth.BigSemiAxis();
   x1 *= norm; y1 *= norm; z1 *= norm;
   x2 *= norm; y2 *= norm; z2 *= norm;

   norm = sqrt(x1*x1 + y1*y1 + z1*z1) * sqrt(x2*x2 + y2*y2 + z2*z2);
   return  (RCore::isEqual(norm, (T)0.)) ? NAN
                           :acos((x1*x2+y1*y2+z1*z2)/norm) * earth.BigSemiAxis();
}





//------------------------------------------------------------------------------
bool _stskmatr(const T &B, const T &L, const T &H,
               T &e11, T &e21, T &e31,
               T &e12, T &e22, T &e32,
               T &e13, T &e23, T &e33,
               T &shX, T &shY, T &shZ, const TEarthDatum &earth)
{
  T  SINL = RCore::_sin(L), COSL = RCore::_cos(L),
     SINB = RCore::_sin(B), COSB = RCore::_cos(B),
     e2   = earth.EllipticityEarthDeriv(),
     tmp  = RCore::_sqrt((T)1.-e2*SINB*SINB);

  if (RCore::isEqual(tmp, (T)0.)) return false;

  T  M    = earth.BigSemiAxis()/tmp;

  shX = (M+H)*COSB*COSL,       //0.
  shY = (M+H)*COSB*SINL,       //(M+H) - e2*M*SINB*SINB
  shZ = (M*((T)1.-e2)+H)*SINB; //-e2*M*SINB*COSB;

  e11 = SINL;
  e21 = COSB*COSL;
  e31 = -SINB*COSL;
  e12 = -COSL;
  e22 = COSB*SINL;
  e32 = -SINB*SINL;
  e13 = 0.;
  e23 = SINB;
  e33 = COSB;

  return true;
}


//------------------------------------------------------------------------------
//Расчёт матрицы перехода от МПСК <=> СТСК
bool mpsk2stskMatr(const T &mB, const T &mL, const T &mH,
                   const T &sB, const T &sL, const T &sH, const T &sQ,
                   T &e11, T &e21, T &e31,
                   T &e12, T &e22, T &e32,
                   T &e13, T &e23, T &e33,
                   T &shX, T &shY, T &shZ, const TEarthDatum &earth)
{
  RMatrix::CPUMatrix M33(3,3),
                     M(1,3);
  if (!_stskmatr(mB, mL, mH,
                 M33.v(0,0), M33.v(1,0), M33.v(2,0),
                 M33.v(0,1), M33.v(1,1), M33.v(2,1),
                 M33.v(0,2), M33.v(1,2), M33.v(2,2),
                 M[0], M[1], M[2], earth))
  {
    return false;
  }

  RMatrix::CPUMatrix S33(3,3),
                     S(1,3);
  if (!_stskmatr(sB, sL, sH,
                 S33.v(0,0), S33.v(1,0), S33.v(2,0),
                 S33.v(0,1), S33.v(1,1), S33.v(2,1),
                 S33.v(0,2), S33.v(1,2), S33.v(2,2),
                 S[0], S[1], S[2], earth))
  {
    return false;
  }

  const T SINQ = RCore::_sin(sQ),
          COSQ = RCore::_cos(sQ);
  T tmp[9] = { -SINQ,  0.,  -COSQ,
              0.,      1.,  0.,
              COSQ,    0.,  -SINQ};

  //матрица доворота старта на север
  RMatrix::CPUMatrix *P33 = RMatrix::CPUMatrix::from2DArray(tmp, 3, 3);
  if (!P33) return false;

  if (!M33.transpose(false)) 
  {
    delete P33;  
    return false;
  }

  
  RMatrix::CPUMatrix _T = M33*S33;
  if (!_T.isValid())
  {
    delete P33;
    return false;
  }
  RMatrix::CPUMatrix _K = S-M;
  if (!_K.isValid())
  {
    delete P33;
    return false;
  }
  RMatrix::CPUMatrix _L = M33*_K;
  if (!_L.isValid())
  {
    delete P33;
    return false;
  }
  RMatrix::CPUMatrix _F = _T*(*P33);
  if (!_F.isValid())
  {
    delete P33;
    return false;
  }
  delete P33;

  e11 = _F.v(0,0);
  e21 = _F.v(1,0);
  e31 = _F.v(2,0);
  e12 = _F.v(0,1);
  e22 = _F.v(1,1);
  e32 = _F.v(2,1);
  e13 = _F.v(0,2);
  e23 = _F.v(1,2);
  e33 = _F.v(2,2);
  shX = _L[0];
  shY = _L[1];
  shZ = _L[2];

  return true;
}



//------------------------------------------------------------------------------
//Преобразование МПСК <=> СТСК (альтернативное преобразование)
bool mpsk2stskA(const T &iX,  const T &iY,  const T &iZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ)
{
  return pgsk2mpsk(iX,iY,iZ,
                   shX,shY,shZ,
                   e11,e21,e31,
                   e12,e22,e32,
                   e13,e23,e33,
                   oX,oY,oZ);
}
bool mpsk2stskA(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ)
{
  return pgsk2mpsk(iX,iY,iZ,iVX,iVY,iVZ,
                   shX,shY,shZ,
                   e11,e21,e31,
                   e12,e22,e32,
                   e13,e23,e33,
                   oX,oY,oZ,oVX,oVY,oVZ);
}

bool stsk2mpskA(const T &iX,  const T &iY,  const T &iZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ)
{
  return mpsk2pgsk(iX,iY,iZ,
                   shX,shY,shZ,
                   e11,e21,e31,
                   e12,e22,e32,
                   e13,e23,e33,
                   oX,oY,oZ);
}
bool stsk2mpskA(const T &iX,  const T &iY,  const T &iZ,
                const T &iVX, const T &iVY, const T &iVZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ)
{
  return mpsk2pgsk(iX,iY,iZ,iVX,iVY,iVZ,
                   shX,shY,shZ,
                   e11,e21,e31,
                   e12,e22,e32,
                   e13,e23,e33,
                   oX,oY,oZ,oVX,oVY,oVZ);
}


//------------------------------------------------------------------------------
//Преобразование координат из АПСК в ОБСК
bool apsk2obsk(const T &x, const T &y, const T &z,
               const T &VX, const T &VY, const T &VZ,
               T &oR, T &oU, T &oV, T &odR, T &odU, T &odV)
{
  if ( RCore::isNan(x) || RCore::isNan(y) || RCore::isNan(z) )
  {
    return false;
  }

  T R = RCore::_sqrt(_POW2(x) + _POW2(y) + _POW2(z));
#ifdef G_ASSERT_SHEME
   assert(R != (T)0.);
#else
   if ( R==(T)0. ) return false;
#endif
  T invR = 1./R;
  T tmp[3] = {x*invR, y*invR, z*invR};

  oR = R;
  oU = tmp[2];
  oV = tmp[1];

  //ПEPECЧET BEKTOPA CKOPOCTИ ИЗ АПСК В OБCK
  if ( (RCore::isNan(VX)==false) && (RCore::isNan(VY)==false) &&
       (RCore::isNan(VZ)==false) )
  {
    odR = (tmp[0]*VX) + (tmp[1]*VY) + (tmp[2]*VZ);
    odU = (VZ - tmp[2]*odR)/R;
    odV = (VY - tmp[1]*odR)/R;
  }
  else
  {
    odR = odU = odV = NAN;
  }
  return true;
}


//------------------------------------------------------------------------------
//Преобразование координат из АПСК в ОБСК
bool obsk2apsk(const T &iR, const T &iU, const T &iV,
               const T &idR, const T &idU, const T &idV,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ)
{
  if ( RCore::isNan(iR) || RCore::isNan(iU) || RCore::isNan(iV) )
  {
    return false;
  }
  T R = 1.0 - _POW2(iU) - _POW2(iV);
#ifdef G_ASSERT_SHEME
  assert(R >= (T)0.);
#else
  if ( R<(T)0. ) return false;
#endif
  R = RCore::_sqrt(R);
#ifdef G_ASSERT_SHEME
  assert(R !=0.);
#else
  if ( R==(T)0. ) return false;
#endif

  oX = R*iR;
  oY = iV*iR;
  oZ = iU*iR;

  //ПEPECЧET BEKTOPA CKOPOCTИ ИЗ OБCK В АПСК
  if ( (RCore::isNan(idR)==false) && (RCore::isNan(idU)==false) &&
       (RCore::isNan(idV)==false) )
  {
    oVX = idR*R - iR/R*(iU*idU + iV*idV);
    oVY = iR*idV + iV*idR;
    oVZ = iR*idU + iU*idR;
  }
  else
  {
    oVX = oVY = oVZ = NAN;
  }
  return true;
}




//------------------------------------------------------------------------------
//Расчёт поправочных значений координат с учётом разницы эпох
bool calcPGSKDeviationsOnEpoch(ruint8 imonth, ruint16 iyear,
                               ruint8 omonth, ruint16 oyear,
                               T &dX, T &dY, T &dZ, const TEarthDatum &earth)
{
  dX = 0.;
  dY = 0.;
  dZ = 0.;
  if ( (imonth>12) || (imonth==0) || (omonth>12) || (omonth==0) ||
       (iyear<2002) || (oyear<2002) ||  ((iyear==oyear) && (imonth==omonth)) )
  {
     return false;
  }

  switch (earth.type())
  {
    case TEarthDatum::ET_PZ90_02:
    case TEarthDatum::ET_PZ90_11:
    case TEarthDatum::ET_GSK2011:
    {
      if ( (iyear<2010) || (oyear<2010) ) return false;
      break;
    }
    case TEarthDatum::ET_PZ90:
    {
      if ( (iyear<2010) || (oyear<2010) ) return false;
      break;
    }
    default: return false;
  }

  T time1 = (imonth<10) ?(T)imonth*0.1 :0.9;
  time1 += (T)iyear;
  T time2 = (omonth<10) ?(T)omonth*0.1 :0.9;
  time2 += (T)oyear;
  T dtime = time2 - time1;
  dX = -0.0212*dtime;
  dY = +0.124 *dtime;
  dZ = +0.0072*dtime;

  return true;
}



void calculateDeviations(const T &iX, const T &iY, const T &iZ,
                         T &oX, T &oY, T &oZ,
                         T &e11, T &e21, T &e31,
                         T &e12, T &e22, T &e32,
                         T &e13, T &e23, T &e33,
                         T &shX, T &shY, T &shZ)
{
  oX = (e11*iX + e21*iY + e31*iZ) + shX;
  oY = (e12*iX + e22*iY + e32*iZ) + shY;
  oZ = (e13*iX + e23*iY + e33*iZ) + shZ;
}

//------------------------------------------------------------------------------
//Формирование матриц коррекции координат для переходов между датумами
rint32 earthDevMatrix(const TEarthDatum &fromEarth, const TEarthDatum &toEarth,
                    T &e11, T &e21, T &e31,
                    T &e12, T &e22, T &e32,
                    T &e13, T &e23, T &e33,
                    T &shX, T &shY, T &shZ)
{
    if (fromEarth==toEarth) return -1;

    switch (fromEarth.type())
    {
      //СК 42
      case TEarthDatum::ET_CK42:
      {
         //СК 42 ==> ПЗ 90
         if (toEarth.type()==TEarthDatum::ET_PZ90)
         {
           e11 = 1.;            e21 = -3.1998e-6;          e31 = 1.6968e-6;
           e12 = 3.1998e-6;     e22 = 1.;                   e32 = 0.;
           e13 = -1.6968e-6;    e23 = 0.;                   e33 = 1.;
           shX  = 25.;          shY = -141.;                shZ = -80.;

           return 1;
         }
         //СК 42 ==> ПЗ 90.02
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           T mult = 1. - 0.22e-6;
           e11 = mult;             e21 = -3.83e-6*mult;   e31 = 1.6968e-6*mult;
           e12 = 3.83e-6*mult;     e22 = mult;            e32 = 0.;
           e13 = -1.6968e-6*mult;  e23 = 0.;              e33 = mult;
           shX  = 23.93;           shY = -141.03;         shZ = -79.98;

           return 1;
         }
         //СК 42 ==> ПЗ 90.11
         else if (toEarth.type()==TEarthDatum::ET_PZ90_11)
         {
           T mult = 1. - 0.228e-6;
           e11 = mult;                       e21 = -3.85043873674005e-6*mult;   e31 = 1.67968547957210e-6*mult;
           e12 = 3.85043873674005e-6*mult;    e22 = mult;                        e32 = -1.11507146655193e-8*mult;
           e13 = -1.67968547957210e-6*mult;  e23 = 1.11507146655193e-8*mult;    e33 = mult;
           shX  = 23.557;                    shY = -140.844;                    shZ = -79.778;

           return 1;
         }
      }
      //ПЗ 90
      case TEarthDatum::ET_PZ90:
      {
         //ПЗ 90 ==> СК 42
         if (toEarth.type()==TEarthDatum::ET_CK42)
         {
           e11 = 1.;            e21 = 3.1998e-6;            e31 = -1.6968e-6;
           e12 = -3.1998e-6;    e22 = 1.;                   e32 = 0.;
           e13 = 1.6968e-6;     e23 = 0.;                   e33 = 1.;
           shX  = -25.;         shY = 141.;                 shZ = 80.;

           return 1;
         }
         //ПЗ 90 ==> ПЗ 90.02
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           T mult = 1. - 0.22e-6;
           e11 = mult;             e21 = -0.6302e-6*mult;   e31 = 0.;
           e12 = 0.6302e-6*mult;   e22 = mult;              e32 = 0.;
           e13 = 0.;               e23 = 0.;                e33 = mult;
           shX  = -1.07;           shY = -0.03;             shZ = 0.02;

           return 1;
         }
         //ПЗ 90 ==> ПЗ 90.11
         else if (toEarth.type()==TEarthDatum::ET_PZ90_11)
         {
           T mult = 1. - 0.228e-6;
           e11 = mult;                       e21 = -6.50668441417108e-7*mult;  e31 = -1.71624043112776e-8*mult;
           e12 = 6.50668441417108e-7*mult;   e22 = mult;                       e32 = -1.11507146655193e-8*mult;
           e13 = 1.71624043112776e-8*mult;   e23 = 1.11507146655193e-8*mult;   e33 = mult;
           shX  = -1.443;                    shY = 0.156;                      shZ = 0.222;

           return 1;
         }
         //ПЗ 90 ==> WGS 84
         else if (toEarth.type()==TEarthDatum::ET_WGS84)
         {
           T mult = 1. - 0.12e-6;
           e11 = mult;               e21 = -0.9696e-6*mult;  e31 = 0.;
           e12 = 0.9696e-6*mult;     e22 = mult;             e32 = 0.;
           e13 = 0.;                 e23 = 0.;               e33 = mult;
           shX  = -1.1;              shY = -0.3;             shZ = -0.9;

           return 1;
         }
      }
      //ПЗ 90.02
      case TEarthDatum::ET_PZ90_02:
      {
         //ПЗ 90.02 ==> СК 42
         if (toEarth.type()==TEarthDatum::ET_CK42)
         {
           T mult = 1. + 0.22e-6;
           e11 = mult;             e21 = 3.83e-6*mult;      e31 = -1.6968e-6*mult;
           e12 = -3.83e-6*mult;    e22 = mult;              e32 = 0.;
           e13 = 1.6968e-6*mult;   e23 = 0.;                e33 = mult;
           shX  = -23.93;          shY = 141.03;            shZ = 79.98;

           return 1;
         }
         //ПЗ 90.02 ==> ПЗ 90
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           T mult = 1. + 0.22e-6;
           e11 = mult;             e21 = 0.6302e-6*mult;    e31 = 0.;
           e12 = -0.6302e-6*mult;  e22 = mult;              e32 = 0.;
           e13 = 0.;               e23 = 0.;                e33 = mult;
           shX = 1.07;           shY = 0.03;             shZ = -0.02;

           return 1;
         }
         //ПЗ 90.02 ==> ПЗ 90.11
         else if (toEarth.type()==TEarthDatum::ET_PZ90_11)
         {
           T mult = 1. - 0.008e-6;
           e11 = mult;                       e21 = -2.04106559747115e-8*mult;  e31 = -1.71624043112776e-8*mult;
           e12 = 2.04106559747115e-8*mult;   e22 = mult;                       e32 = -1.11507146655193e-8*mult;
           e13 = 1.71624043112776e-8*mult;   e23 = 1.11507146655193e-8*mult;   e33 = mult;
           shX  = -0.373;                    shY = 0.186;                      shZ = 0.202;

           return 1;
         }
         //ПЗ 90.02 ==> WGS 84
         else if (toEarth.type()==TEarthDatum::ET_WGS84)
         {
           e11 = 1.;      e21 = 0.;    e31 = 0.;
           e12 = 0.;      e22 = 1.;    e32 = 0.;
           e13 = 0.;      e23 = 0.;    e33 = 1.;
           shX  = -0.36;  shY = 0.08;  shZ = 0.18;

           return 1;
         }
         //ПЗ 90.02 ==> ГСК 2011
         else if (toEarth.type()==TEarthDatum::ET_GSK2011)
         {
           T mult = 1. + 0.0086e-6;
           e11 = mult;                       e21 = -2.06676072256995e-8*mult;  e31 = -1.72545189106884e-8*mult;
           e12 = 2.06676072256995e-8*mult;    e22 = mult;                       e32 = -8.42606177768373e-9*mult;
           e13 = 1.72545189106884e-8*mult;   e23 = 8.42606177768373e-9*mult;   e33 = mult;
           shX  = -0.373;                    shY = 0.172;                      shZ = 0.21;

           return 1;
         }
      }
      //ПЗ 90.11
      case TEarthDatum::ET_PZ90_11:
      {
         //ПЗ 90.11 ==> СК 42
         if (toEarth.type()==TEarthDatum::ET_CK42)
         {
           T mult = 1. + 0.228e-6;
           e11 = mult;                       e21 = 3.85043873674005e-6*mult;    e31 = -1.67968547957210e-6*mult;
           e12 = -3.85043873674005e-6*mult;  e22 = mult;                        e32 = 1.11507146655193e-8*mult;
           e13 = 1.67968547957210e-6*mult;   e23 = -1.11507146655193e-8*mult;   e33 = mult;
           shX = -23.557;                    shY = 140.844;                     shZ = 79.778;

           return 1;
         }
         //ПЗ 90.11 ==> ПЗ 90.02
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           T mult = 1. + 0.008e-6;
           e11 = mult;                       e21 = 2.04106559747115e-8*mult;    e31 = 1.71624043112776e-8*mult;
           e12 = -2.04106559747115e-8*mult;  e22 = mult;                        e32 = 1.11507146655193e-8*mult;
           e13 = -1.71624043112776e-8*mult;  e23 = -1.11507146655193e-8*mult;   e33 = mult;
           shX = 0.373;                      shY = -0.186;                      shZ = -0.202;

           return 1;
         }
         //ПЗ 90.11 ==> ПЗ 90
         else if (toEarth.type()==TEarthDatum::ET_PZ90)
         {
           T mult = 1. + 0.228e-6;
           e11 = mult;                       e21 = 6.50668441417108e-7*mult;  e31 = 1.71624043112776e-8*mult;
           e12 = -6.50668441417108e-7*mult;  e22 = mult;                      e32 = 1.11507146655193e-8*mult;
           e13 = -1.71624043112776e-8*mult;  e23 = -1.11507146655193e-8*mult; e33 = mult;
           shX  = 1.443;                     shY = -0.156;                    shZ = -0.222;

           return 1;
         }
         //ПЗ 90.11 ==> WGS 84
         else if (toEarth.type()==TEarthDatum::ET_WGS84)
         {
           T mult = 1. + 0.008e-6;
           e11 = mult;                       e21 = 2.04106559747115e-8*mult;    e31 = 1.71624043112776e-8*mult;
           e12 = -2.04106559747115e-8*mult;  e22 = mult;                        e32 = 1.11507146655193e-8*mult;
           e13 = -1.71624043112776e-8*mult;  e23 = -1.11507146655193e-8*mult;   e33 = mult;
           shX = 0.013;                      shY = -0.106;                      shZ = -0.022;

           return 1;
         }
         //ПЗ 90.11 ==> ГСК 2011
         else if (toEarth.type()==TEarthDatum::ET_GSK2011)
         {
           T mult = 1. + 0.0006e-6;
           e11 = mult;                       e21 = -2.56951250988054e-10*mult; e31 = -9.21145994108119e-11*mult;
           e12 = 2.56951250988054e-10*mult;  e22 = mult;                       e32 = 2.72465288783559e-9*mult;
           e13 = 9.21145994108119e-11*mult;  e23 = -2.72465288783559e-9*mult;  e33 = mult;
           shX  = 0.;                        shY = -0.014;                     shZ = 0.008;

           return 1;
         }
      }
      //ПЗ WGS 84
      case TEarthDatum::ET_WGS84:
      {
         //WGS 84 ==> ПЗ 90
         if (toEarth.type()==TEarthDatum::ET_PZ90)
         {
           T mult = 1. + 0.12e-6;
           e11 = mult;              e21 = -0.9696e-6*mult;  e31 = 0.;
           e12 = 0.9696e-6*mult;    e22 = mult;             e32 = 0.;
           e13 = 0.;                e23 = 0.;               e33 = mult;
           shX  = 1.1;              shY = 0.3;              shZ = 0.9;

           return 1;
         }
         //WGS 84 ==> ПЗ 90.02
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           e11 = 1.;              e21 = 0.;                 e31 = 0.;
           e12 = 0.;              e22 = 1.;                 e32 = 0.;
           e13 = 0.;              e23 = 0.;                 e33 = 1.;
           shX = 0.36;            shY = -0.08;              shZ = -0.18;

           return 1;
         }
         //WGS 84 ==> ПЗ 90.11
         else if (toEarth.type()==TEarthDatum::ET_PZ90_11)
         {
           T mult = 1. - 0.008e-6;
           e11 = mult;                       e21 = -2.04106559747115e-8*mult;   e31 = -1.71624043112776e-8*mult;
           e12 = 2.04106559747115e-8*mult;   e22 = mult;                        e32 = -1.11507146655193e-8*mult;
           e13 = 1.71624043112776e-8*mult;   e23 = 1.11507146655193e-8*mult;    e33 = mult;
           shX = -0.013;                     shY = 0.106;                       shZ = 0.022;

           return 1;
         }
         //WGS 84 ==> ГСК 2011
         else if (toEarth.type()==TEarthDatum::ET_GSK2011)
         {
           T mult = 1. + 0.0086e-6;
           e11 = mult;                      e21 = -2.06676072256995e-8*mult; e31 = -1.72545189106884e-8*mult;
           e12 = 2.06676072256995e-8*mult;  e22 = mult;                      e32 = -8.42606177768373e-9*mult;
           e13 = 1.72545189106884e-8*mult;  e23 = 8.42606177768373e-9*mult;  e33 = mult;
           shX  = -0.013;                   shY = 0.092;                     shZ = 0.003;

           return 1;
         }
      }
      //ПЗ ГСК 2011
      case TEarthDatum::ET_GSK2011:
      {
         //ГСК 2011 ==> WGS 84
         if (toEarth.type()==TEarthDatum::ET_WGS84)
         {
           T mult = 1. - 0.0086e-6;
           e11 = mult;                       e21 = 2.06676072256995e-8*mult;  e31 = 1.72545189106884e-8*mult;
           e12 = -2.06676072256995e-8*mult;  e22 = mult;                      e32 = 8.42606177768373e-9*mult;
           e13 = -1.72545189106884e-8*mult;  e23 = -8.42606177768373e-9*mult; e33 = mult;
           shX  = 0.013;                     shY = -0.092;                    shZ = -0.03;

           return 1;
         }
         //ГСК 2011 ==> ПЗ 90.02
         else if (toEarth.type()==TEarthDatum::ET_PZ90_02)
         {
           T mult = 1. - 0.0086e-6;
           e11 = mult;                       e21 = 2.06676072256995e-8*mult;  e31 = 1.72545189106884e-8*mult;
           e12 = -2.06676072256995e-8*mult;  e22 = mult;                      e32 = 8.42606177768373e-9*mult;
           e13 = -1.72545189106884e-8*mult;  e23 = -8.42606177768373e-9*mult; e33 = mult;
           shX  = 0.373;                     shY = -0.172;                    shZ = -0.021;

           return 1;
         }
         //ГСК 2011 ==> ПЗ 90.11
         else if (toEarth.type()==TEarthDatum::ET_PZ90_11)
         {
           T mult = 1. - 0.006e-6;
           e11 = mult;                       e21 = 2.56951250988054e-10*mult;   e31 = 9.21145994108119e-11*mult;
           e12 = -2.56951250988054e-10*mult; e22 = mult;                        e32 = -2.72465288783559e-09*mult;
           e13 = -9.21145994108119e-11*mult; e23 = 2.72465288783559e-09*mult;   e33 = mult;
           shX = 0.;                         shY = 0.014;                       shZ = -0.008;

           return 1;
         }
      }

      default: return 0;
    }
}




} //namespace RGeoMath

} //namespace RMath

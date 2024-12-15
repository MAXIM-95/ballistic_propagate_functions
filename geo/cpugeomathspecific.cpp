#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cpumath.h"
#include "cpumatrix.h"
#include "cpugeomath.h"
#include "cpugeomathspecific.h"

namespace RMath {



namespace RGeoMath {



bool stskmpskInternal(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               size_t ichannel_number, T *TR, T *SH, bool two)
{
  if ( RCore::isNan(iX) || RCore::isNan(iY) || RCore::isNan(iZ) )
  {
    return false;
  }

  static const size_t first_channel_index = 1;
#ifdef G_ASSERT_SHEME
  assert(ichannel_number>=first_channel_index);
#else
  if (ichannel_number<first_channel_index) return false;
#endif

  ruint32 sft = ichannel_number - first_channel_index,
          id1 = sft*3,
          id2 = sft*9;

  return (two) ?pgsk2mpsk(iX, iY, iZ, iVX, iVY, iVZ,
                          SH[id1],   SH[id1+1], SH[id1+2],
                          TR[id2],   TR[id2+1], TR[id2+2],
                          TR[id2+3], TR[id2+4], TR[id2+5],
                          TR[id2+6], TR[id2+7], TR[id2+8],
                          oX, oY, oZ, oVX, oVY, oVZ)
               :mpsk2pgsk(iX, iY, iZ, iVX, iVY, iVZ,
                          SH[id1],   SH[id1+1], SH[id1+2],
                          TR[id2],   TR[id2+1], TR[id2+2],
                          TR[id2+3], TR[id2+4], TR[id2+5],
                          TR[id2+6], TR[id2+7], TR[id2+8],
                          oX, oY, oZ, oVX, oVY, oVZ);
}


bool apskmpskInternal(const T &iX,  const T &iY,  const T &iZ,
                      const T &iVX, const T &iVY, const T &iVZ,
                      T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
                      size_t nEdge, T *TR, T *SH, bool two, ruint32 edge_count)
{
  if ( RCore::isNan(iX) || RCore::isNan(iY) || RCore::isNan(iZ) )
  {
    return false;
  }

  static const size_t first_edge_index = 0;
#ifdef G_ASSERT_SHEME
  assert(nEdge<edge_count);
#else
  if (nEdge >= edge_count) return false;
#endif

  ruint32 sft = nEdge - first_edge_index,
          id1 = sft*3,
          id2 = sft*9;

  return (two) ?mpsk2pgsk(iX, iY, iZ, iVX, iVY, iVZ,
                          SH[id1],   SH[id1+1], SH[id1+2],
                          TR[id2],   TR[id2+3], TR[id2+6],
                          TR[id2+1], TR[id2+4], TR[id2+7],
                          TR[id2+2], TR[id2+5], TR[id2+8],
                          oX, oY, oZ, oVX, oVY, oVZ)
               :pgsk2mpsk(iX, iY, iZ, iVX, iVY, iVZ,
                          SH[id1],   SH[id1+1], SH[id1+2],
                          TR[id2],   TR[id2+3], TR[id2+6],
                          TR[id2+1], TR[id2+4], TR[id2+7],
                          TR[id2+2], TR[id2+5], TR[id2+8],
                          oX, oY, oZ, oVX, oVY, oVZ);
}



namespace GIP10 {


//------------------------------------------------------------------------------
//МАТРИЦЫ ПЕРЕХОДА ДЛЯ СТАРТОВЫХ СК                              |
T TRANSFORM[2*9] =
{

};




//------------------------------------------------------------------------------
// Вектора сдвига из СТСК                            |
T  SHIFT[2*3] =
{

};





bool stsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               size_t ichannel_number)
{
  return stskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, ichannel_number,
                          TRANSFORM, SHIFT, true);
}



bool mpsk2stsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               size_t ichannel_number)
{
  return stskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, ichannel_number,
                          TRANSFORM, SHIFT, false);
}






const T ca = RCore::_cos(a),
        sa = RCore::_sin(a),
        ce = RCore::_cos(e),
        se = RCore::_sin(e),
        cb = RCore::_cos(b),
        sb = RCore::_sin(b);

T	MRLSSHIFT[3]     = { -ac*sa - bc*ca,     -hc,    ac*ca - bc*sa };
T	MRLSTRANSFORM[9] = { -se*sa,              ce,    se*ca,
                              cb*ce*sa - sb*ca,   cb*se, -cb*ce*ca - sb*sa,
                             -sb*ce*sa - cb*ca,  -sb*se,  sb*ce*ca - cb*sa };

//------------------------------------------------------------------------------
//Пересчет из АПСК в МПСК                                                      |
bool apsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ, ruint32 nEdge)
{
  return apskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, nEdge,
                          MRLSTRANSFORM, MRLSSHIFT, true, 1);
}


//------------------------------------------------------------------------------
//Пересчет из МПСК в АПСК                                                      |
bool mpsk2apsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ, ruint32 nEdge)
{
  return apskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, nEdge,
                          MRLSTRANSFORM, MRLSSHIFT, false, 1);
}





}





namespace SOFRINO {


//------------------------------------------------------------------------------
//МАТРИЦЫ ПЕРЕХОДА ДЛЯ СТАРТОВЫХ СК В/ИЗ МПСК МРЛС                             |
T TRANSFORM[68*9] = {


};




//------------------------------------------------------------------------------
// Вектора сдвига из СТСК в МПК (начало МПК в СТСК)                            |
T  SHIFT[68*3] = {


};



//------------------------------------------------------------------------------
// МАТРИЦЫ ПЕРЕСЧЕТА ДЛЯ СИСТЕМ МПК И ОБСК                                     |
T MLPER[36] =
{


};


T  MLBPER[12] =
{


};

// параметры для выбора номера
real64 APER = 0;
real64 BPER = 0.;







//------------------------------------------------------------------------------
//Пересчет из Стартовой СК в МПСК                                              |
//вх.данные:  [iX,iY,iZ] - координаты вектора положения в Стартовой СК         |
//            [iVX,iVY,iVZ] - координаты вектора скорости в Стартовой СК       |
//вых.данные: [oX,oY,oZ,oVX,oVY,oVZ] МПСК координаты [м] и их производные[м/с] |
bool stsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               size_t ichannel_number)
{
  return stskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, ichannel_number,
                          TRANSFORM, SHIFT, true);
}


//------------------------------------------------------------------------------
//Пересчет из МПСК в Стартовую СК                                              |
//вх.данные:  [iX,iY,iZ] - координаты вектора положения в МПСК                 |
//            [iVX,iVY,iVZ] - координаты вектора скорости в МПСК                |
//вых.данные: [oX,oY,oZ,oVX,oVY,oVZ] координаты [м] и их производные[м/с]      |
//                                   в Стартовой СК                            |
bool mpsk2stsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               size_t ichannel_number)
{
  return stskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, ichannel_number,
                          TRANSFORM, SHIFT, false);
}




//------------------------------------------------------------------------------
//Пересчет из ОБСК в МПСК                                                      |
bool obsk2mpsk(const T &iR,  const T &iU,  const T &iV,
               const T &idR, const T &idU, const T &idV,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ,
               rint32 i_NK, rint32 i_startNK, rint32 nEdge)
{
  if ( RCore::isNan(iR) || RCore::isNan(iU) || RCore::isNan(iV) )
  {
    return false;
  }

  T R1  =  (T)1. - iU*iU - iV*iV;
#ifdef G_ASSERT_SHEME
  assert(R1>=(T)0.);
#else
  if (R1 < (T)0.) return false;
#endif
  R1 = RCore::_sqrt(R1);

  T XGL = iR*R1,
    YGL = iR*iV,
    ZGL = iR*iU;
  //матрица пересчета из МПК в ОБСК (из МПК в АСПК) в зависимости от грани МРЛС
  T ML[9],
  //вектор пересчета из МПК в ОБСК (из МПК в АСПК) в зависимости от грани МРЛС
    MLV[3];

  if ( (i_startNK==1 && i_NK>68) ||
       (i_startNK==33 &&  ((i_NK<33) || i_NK>(33+68-1))) ) return false;

  //Выбор матрицы пересчета
  //запись на ЭЛЬ-76 (для справки):   ML:=MLПEP[(L-1)*9:9];
  for (rint32 i=0; i<9; i++)
    ML[i] = MLPER[9*(nEdge-1)+i];

  //Выбор вектора пересчета
  //запись на ЭЛЬ-76 (для справки): MLV:=MLBПEP[(L-1)*3:3];
  for (rint32 i=0; i<3; i++)
    MLV[i] = MLBPER[3*(nEdge-1)+i];

  //ПEPECЧET BEKTOPA ПOЛOЖEHИЯ ИЗ OБCK B MПCK
  XGL -= MLV[0];
  YGL -= MLV[1];
  ZGL -= MLV[2];
  oX = ML[0]*XGL + ML[3]*YGL + ML[6]*ZGL;
  oY = ML[1]*XGL + ML[4]*YGL + ML[7]*ZGL;
  oZ = ML[2]*XGL + ML[5]*YGL + ML[8]*ZGL;

  //ПEPECЧET BEKTOPA CKOPOCTИ ИЗ OБCK B MПCK
  if ( (RCore::isNan(idR)==false) && (RCore::isNan(idU)==false) &&
       (RCore::isNan(idV)==false) )
  {
    XGL = idR*R1 - iR*(iU*idU + iV*idV)/R1;
    YGL = iR*idV + idR*iV;
    ZGL = iR*idU + idR*iU;
    oVX = ML[0]*XGL + ML[3]*YGL + ML[6]*ZGL;
    oVY = ML[1]*XGL + ML[4]*YGL + ML[7]*ZGL;
    oVZ = ML[2]*XGL + ML[5]*YGL + ML[8]*ZGL;
  }
  else
  {
    oVX = oVY = oVZ = NAN;
  }

  return true;
}



//------------------------------------------------------------------------------
//Пересчет из МПСК в ОБСК                                                      |
bool mpsk2obsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oR, T &oU, T &oV, T &odR, T &odU,
               T &odV, rint32 i_NK, rint32 i_startNK, rint32 &nEdge)
{
  if ( RCore::isNan(iX) || RCore::isNan(iY) || RCore::isNan(iZ) )
  {
    return false;
  }

  T XGLT = 0.f, YGLT = 0.f, ZGLT = 0.f,
    XGL  = 0.f, YGL  = 0.f, ZGL  = 0.f,
    RYA  = APER*iX + iZ,
    RYA1 = -BPER*iX + iZ;
  //матрица пересчета из МПК в ОБСК (из МПК в АСПК) в зависимости от грани МРЛС
  T ML[9],
  //вектор пересчета из МПК в ОБСК (из МПК в АСПК) в зависимости от грани МРЛС
    MLV[3];

  if ( (i_startNK==1 && i_NK>68) ||
       (i_startNK==33 && ((i_NK<33) || i_NK>(33+68-1))) ) return false;

  rint32 iedge;
  //Выбор номера грани МРЛС
  if (RYA<=0.f && RYA1>0.f) iedge = 3;
  else if (RYA<0.f && RYA1<=0.f) iedge = 4;
  else if (RYA>=0.f && RYA1<0.f) iedge = 1;
  else iedge = 2;
  if (nEdge<0) iedge = _MFASTABS(nEdge);
  else nEdge = iedge;

  //Выбор матрицы пересчета
  for (rint32 i=0; i<9; i++)
    ML[i] = MLPER[9*(iedge-1) + i];

  //Выбор вектора пересчета

  for (rint32 i=0; i<3; i++)
    MLV[i] = MLBPER[3*(iedge-1) + i];

  XGL = ML[0]*iX + ML[1]*iY + ML[2]*iZ + MLV[0];
  YGL = ML[3]*iX + ML[4]*iY + ML[5]*iZ + MLV[1];
  ZGL = ML[6]*iX + ML[7]*iY + ML[8]*iZ + MLV[2];

  oR = RCore::_sqrt(XGL*XGL + YGL*YGL + ZGL*ZGL);
  if (RCore::isEqual(oR, (T)0.)) return false;

  //ПEPECЧET BEKTOPA ПOЛOЖEHИЯ ИЗ MПK B OБCK
  oV = YGL/oR;
  oU = ZGL/oR;

  //ПEPECЧET BEKTOPA CKOPOCTИ ИЗ OБCK B MПK
  if ( (RCore::isNan(iVX)==false) && (RCore::isNan(iVY)==false) &&
       (RCore::isNan(iVZ)==false) )
  {
    XGLT = ML[0]*iVX + ML[1]*iVY + ML[2]*iVZ;
    YGLT = ML[3]*iVX + ML[4]*iVY + ML[5]*iVZ;
    ZGLT = ML[6]*iVX + ML[7]*iVY + ML[8]*iVZ;
    odR  = (XGL*XGLT + YGL*YGLT + ZGL*ZGLT)/oR;
    odV  = (oR*YGLT - YGL*odR) / (oR*oR);
    odU  = (oR*ZGLT - ZGL*odR) / (oR*oR);
  }
  else
  {
    odR = odV = odU = NAN;
  }

  return true;
}





const ruint32 edge_cnt = 4;

//Азимуты нормалей к граням
T a[edge_cnt] =
{

};
//Углы между горизонтальной плоскостью и антеннами граней
T e[edge_cnt] =
{

};
//Поворот приемной ФАР4
T b[edge_cnt] =
{

};


const T ac = 0, // Полурасстояние между гранями
        bc = 0,// Сдвиг
        hc = 0;// Полувысота


T	MRLSSHIFT[3*edge_cnt] =
{
   -ac*RCore::_sin(a[0]) - bc*RCore::_cos(a[0]),   -hc,
    ac*RCore::_cos(a[0]) - bc*RCore::_sin(a[0]),

   -ac*RCore::_sin(a[1]) - bc*RCore::_cos(a[1]),   -hc,
    ac*RCore::_cos(a[1]) - bc*RCore::_sin(a[1]),

   -ac*RCore::_sin(a[2]) - bc*RCore::_cos(a[2]),   -hc,
    ac*RCore::_cos(a[2]) - bc*RCore::_sin(a[2]),

   -ac*RCore::_sin(a[3]) - bc*RCore::_cos(a[3]),   -hc,
    ac*RCore::_cos(a[3]) - bc*RCore::_sin(a[3])
};

T	MRLSTRANSFORM[9*edge_cnt] =
{
  -RCore::_sin(e[0])*RCore::_sin(a[0]),
   RCore::_cos(e[0]),
   RCore::_sin(e[0])*RCore::_cos(a[0]),
   RCore::_cos(b[0])*RCore::_cos(e[0])*RCore::_sin(a[0]) - RCore::_sin(b[0])*RCore::_cos(a[0]),
   RCore::_cos(b[0])*RCore::_sin(e[0]),
  -RCore::_cos(b[0])*RCore::_cos(e[0])*RCore::_cos(a[0]) - RCore::_sin(b[0])*RCore::_sin(a[0]),
  -RCore::_sin(b[0])*RCore::_cos(e[0])*RCore::_sin(a[0]) - RCore::_cos(b[0])*RCore::_cos(a[0]),
  -RCore::_sin(b[0])*RCore::_sin(e[0]),
   RCore::_sin(b[0])*RCore::_cos(e[0])*RCore::_cos(a[0]) - RCore::_cos(b[0])*RCore::_sin(a[0]),

  -RCore::_sin(e[1])*RCore::_sin(a[1]),
   RCore::_cos(e[1]),
   RCore::_sin(e[1])*RCore::_cos(a[1]),
   RCore::_cos(b[1])*RCore::_cos(e[1])*RCore::_sin(a[1]) - RCore::_sin(b[1])*RCore::_cos(a[1]),
   RCore::_cos(b[1])*RCore::_sin(e[1]),
  -RCore::_cos(b[1])*RCore::_cos(e[1])*RCore::_cos(a[1]) - RCore::_sin(b[1])*RCore::_sin(a[1]),
  -RCore::_sin(b[1])*RCore::_cos(e[1])*RCore::_sin(a[1]) - RCore::_cos(b[1])*RCore::_cos(a[1]),
  -RCore::_sin(b[1])*RCore::_sin(e[1]),
   RCore::_sin(b[1])*RCore::_cos(e[1])*RCore::_cos(a[1]) - RCore::_cos(b[1])*RCore::_sin(a[1]),

  -RCore::_sin(e[2])*RCore::_sin(a[2]),
   RCore::_cos(e[2]),
   RCore::_sin(e[2])*RCore::_cos(a[2]),
   RCore::_cos(b[2])*RCore::_cos(e[2])*RCore::_sin(a[2]) - RCore::_sin(b[2])*RCore::_cos(a[2]),
   RCore::_cos(b[2])*RCore::_sin(e[2]),
  -RCore::_cos(b[2])*RCore::_cos(e[2])*RCore::_cos(a[2]) - RCore::_sin(b[2])*RCore::_sin(a[2]),
  -RCore::_sin(b[2])*RCore::_cos(e[2])*RCore::_sin(a[2]) - RCore::_cos(b[2])*RCore::_cos(a[2]),
  -RCore::_sin(b[2])*RCore::_sin(e[2]),
   RCore::_sin(b[2])*RCore::_cos(e[2])*RCore::_cos(a[2]) - RCore::_cos(b[2])*RCore::_sin(a[2]),

  -RCore::_sin(e[3])*RCore::_sin(a[3]),
   RCore::_cos(e[3]),
   RCore::_sin(e[3])*RCore::_cos(a[3]),
   RCore::_cos(b[3])*RCore::_cos(e[3])*RCore::_sin(a[3]) - RCore::_sin(b[3])*RCore::_cos(a[3]),
   RCore::_cos(b[3])*RCore::_sin(e[3]),
  -RCore::_cos(b[3])*RCore::_cos(e[3])*RCore::_cos(a[3]) - RCore::_sin(b[3])*RCore::_sin(a[3]),
  -RCore::_sin(b[3])*RCore::_cos(e[3])*RCore::_sin(a[3]) - RCore::_cos(b[3])*RCore::_cos(a[3]),
  -RCore::_sin(b[3])*RCore::_sin(e[3]),
   RCore::_sin(b[3])*RCore::_cos(e[3])*RCore::_cos(a[3]) - RCore::_cos(b[3])*RCore::_sin(a[3])
};

//------------------------------------------------------------------------------
//Пересчет из АПСК в МПСК                                                      |
bool apsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ, ruint32 nEdge)
{
  return apskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, nEdge,
                          MRLSTRANSFORM, MRLSSHIFT, true, edge_cnt);
}


//------------------------------------------------------------------------------
//Пересчет из МПСК в АПСК                                                      |
bool mpsk2apsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ, ruint32 nEdge)
{
  return apskmpskInternal(iX, iY, iZ, iVX, iVY, iVZ,
                          oX, oY, oZ, oVX, oVY, oVZ, nEdge,
                          MRLSTRANSFORM, MRLSSHIFT, false, edge_cnt);
}



}


















} //namespace RGeoMath

} //namespace RMath

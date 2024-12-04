#ifdef G_ASSERT_SHEME
  #include <assert.h>
#endif

#include "cputargeterkoshi.h"
#include "cpuMSIGMballistics.h"

using namespace RModels::RAtmosphere;

namespace RMath {


namespace RBallistics {




//------------------------------------------------------------------------------
//Коррекция требуемых скоростей для попадания в цель с заданной точностью
bool CPUTargeterKoshi::adjustConditions(T &t, const T &x, const T &y, const T &z,
                                        T &vx, T &vy, T &vz,
                                        const T &tx, const T &ty, const T &tz,
                                        const T &bk1, const T &bk2,
                                        const T& idt, T& idelt, ruint32 maxiter,
                                        const T &atmorange)
{
  if ( (!_prop) || RCore::isNan(x) || RCore::isNan(y) || RCore::isNan(z) ||
       RCore::isNan(t) || RCore::isNan(bk1) || RCore::isNan(bk2) ||
       RCore::isNan(vx) || RCore::isNan(vy) || RCore::isNan(vz) ||
       RCore::isNan(tx) || RCore::isNan(ty) || RCore::isNan(tz) ||
       RCore::isNan(idt) || RCore::isNan(idelt) || RCore::isEqual(idt, 0.) ||
       (idelt<=0.) )
  {
   return false;
  }


  T  dVx, dVy, dVz, addx, addy, addz;
  T  iVx  = vx,
     iVy  = vy,
     iVz  = vz;
  T  rrt  = t,
     rVx  = vx,
     rVy  = vy,
     rVz  = vz,
     rdif = NAN;

  CPUMatrix   results(4,3),
              delt(3,3),
              deltx(3,3),
              delty(3,3),
              deltz(3,3),
              targ(1,3);
  targ[0] = tx;
  targ[1] = ty;
  targ[2] = tz;

  T rt = NAN, rx, ry, rz, vrx, vry, vrz;
  T dif = NAN;

  ruint32 iterations_cnt = 0;
  bool flag = true;

  //массив параметров оптимальных НУ
  double lmin[5];
  lmin[0] = (double)NAN;
  lmin[1] = NAN;
  lmin[2] = NAN;
  lmin[3] = NAN;
  lmin[4] = NAN;

  T heighRange = _prop->heightInBaseSystem(tx, ty, tz);


  while (flag && (iterations_cnt < maxiter) )
  {
     //определение реакции системы на единичное воздействие
     for (ruint32 i=0; i<4; i++)
     {
        rx = x;
        ry = y;
        rz = z;
        switch (i)
        {
          case 1:
          {
            vrx = iVx+1;
            vry = iVy;
            vrz = iVz;     //dx
            break;
          }
          case 2:
          {
            vrx = iVx;
            vry = iVy+1;
            vrz = iVz;     //dy
            break;
          }
          case 3:
          {
            vrx = iVx;
            vry = iVy;
            vrz = iVz+1;   //dz
            break;
          }
          default:
          {
            vrx = iVx;
            vry = iVy;
            vrz = iVz;     //0
          }
        }


       rt = 0.;
       dVx = dVy = dVz = 0.;
       CPUBasePropagator::PropagateResult pres =
                    _prop->propagateOnHeight(rt, idt, rx, ry, rz, vrx, vry, vrz,
                                             dVx, dVy, dVz, bk1, bk2, heighRange,
                                             CPUBasePropagator::FD_DESCENDING,
                                             3600., atmorange);

       //ошибка при прогнозе
       if (pres>CPUBasePropagator::PR_OK_WITH_LIMITATIONS)
       {
         flag = false;
         break;
       }

       results.v(i,0) = rx;
       results.v(i,1) = ry;
       results.v(i,2) = rz;

       //поиск оптимального управления
       if (i==0)
       {
          dif = RCore::_sqrt(_POW2(tx-rx) + _POW2(ty-ry) + _POW2(tz-rz));
          if (RCore::isNan(dif)) return false;

          if (RCore::isNan(lmin[4]) || (dif<lmin[4]) )
          {
            lmin[0] = rt;
            lmin[1] = iVx;
            lmin[2] = iVy;
            lmin[3] = iVz;
            lmin[4] = dif;
            //сохраняем исходное отклонение для начальных условий до коррекции
            if (RCore::isNan(rdif)) rdif = dif;
          }
       }

       if ( (i==0) && (dif<=idelt) )
       {
         flag = false;
         break;
       }
       //формируем главный определитель
       else if (i>0)
       {
         for (ruint32 k=0; k<3; k++)
         {
           delt.v(i-1, k) = results.v(i, k) - results.v(0, k);
         }
       }
     }

     if (flag)//формируем матрицы откликов если это необходимо
     {
       deltx = delt; delty = delt; deltz = delt;
       CPUMatrix *rs = results.column(0);
       CPUMatrix dop = targ - (*rs);
       delete rs;

       for (ruint32 k=0; k<3; k++)
       {
          deltx.v(0,k) = dop[k];
          delty.v(1,k) = dop[k];
          deltz.v(2,k) = dop[k];
       }

       addx = deltx.det()/delt.det();
       addy = delty.det()/delt.det();
       addz = deltz.det()/delt.det();
       iVx += addx; iVy += addy; iVz += addz;

       iterations_cnt++;
     }
  }

  bool havent_currency = (RCore::isNan(lmin[4])==true) || (dif<lmin[4]);
  t    += (havent_currency) ?rt  :lmin[0];
  vx    = (havent_currency) ?iVx :lmin[1];
  vy    = (havent_currency) ?iVy :lmin[2];
  vz    = (havent_currency) ?iVz :lmin[3];
  T odelt = (havent_currency) ?dif :lmin[4];

  //если пристрелка только ухудшила результаты
  if (rdif<odelt)
  {
    t  = rrt;
    vx = rVx;
    vy = rVy;
    vz = rVz;
    odelt = rdif;
  }

  flag = odelt<=idelt;
  idelt = odelt;
  return flag;
}




} //namespace RBallistics

} //namespace RMath

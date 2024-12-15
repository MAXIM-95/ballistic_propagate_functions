#pragma once

#include "models/datums/cpugeodatums.h"

using namespace RModels::RDatums;

namespace RMath {


namespace RGeoMath {



//------------------------------------------------------------------------------
//                ФУНКЦИИ ПРЕОБРАЗОВАНИЯ СИСТЕМ КООРДИНАТ                      |
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
/*!
 * \brief Преобразование Геодезических координат в Геоцентрические
 * \param B геодезическая широта, [рад]
 * \param L геодезическая долгота, [рад]
 * \param H геодезическая высота, [м]
 * \param Bgc Геоцентрич. широта, [рад]
 * \param Lgc Геоцентрич. долгота, [рад]
 * \param Hgc Геоцентрич. высота, [м]
 * \param earth Параметры Земли
 */
bool geodez2geocentrical(const T &B,  const T &L,  const T &H,
                         T &Bgc, T &Lgc, T &Hgc, const TEarthDatum &earth=CK42);
/*!
 * \brief Преобразование Геодезических координат в Геоцентрические
 * \param B геодезическая широта, [рад]
 * \param L геодезическая долгота, [рад]
 * \param H геодезическая высота, [м]
 * \param dB Геодез. производная широты, [рад/с]
 * \param dL Геодез. производная долготы, [рад/с]
 * \param dH Геодез. производная высоты, [м/с]
 * \param Bgc Геоцентрич. широта, [рад]
 * \param Lgc Геоцентрич. долгота, [рад]
 * \param Hgc Геоцентрич. высота, [м]
 * \param dBgc Геоц. произв. широты, [рад/c]
 * \param dLgc Геоц. произв. долготы, [рад/c]
 * \param dHgc Геоц. произв. высоты, [м/c]
 * \param earth Параметры Земли
 */
bool geodez2geocentrical(const T &B,  const T &L,  const T &H,
                         const T &dB, const T &dL, const T &dH,
                         T &Bgc, T &Lgc, T &Hgc, T &dBgc, T &dLgc, T &dHgc,
                         const TEarthDatum &earth=CK42);
//------------------------------------------------------------------------------
/*!
 * \brief Преобразование Геоцентрических координат в Геодезические
 * \param Bgc Геоцентрич. широта, [рад]
 * \param Lgc Геоцентрич. долгота, [рад]
 * \param Hgc Геоцентрич. высота, [м]
 * \param B геодезическая широта, [рад]
 * \param L геодезическая долгота, [рад]
 * \param H высота широта, [м]
 * \param earth Параметры Земли
 */
bool geocentrical2geodez(const T &Bgc,  const T &Lgc,  const T &Hgc,
                         T &B, T &L, T &H, const TEarthDatum &earth=CK42);
/*!
 * \brief Преобразование Геоцентрических координат в Геодезические
 * \param Bgc геодезическая широта, [рад]
 * \param Lgc геодезическая долгота, [рад]
 * \param Hgc геодезическая высота, [м]
 * \param dBgc Геодез. производная широты, [рад/с]
 * \param dLgc Геодез. производная долготы, [рад/с]
 * \param dHgc Геодез. производная высоты, [м/с]
 * \param B Геоцентрич. широта, [рад]
 * \param L Геоцентрич. долгота, [рад]
 * \param H Геоцентрич. высота, [м]
 * \param dB Геоц. произв. широты, [рад/c]
 * \param dL Геоц. произв. долготы, [рад/c]
 * \param dH Геоц. произв. высоты, [м/c]
 * \param earth Параметры Земли
 */
bool geocentrical2geodez(const T &Bgc,  const T &Lgc,  const T &Hgc,
                         const T &dBgc,  const T &dLgc,  const T &dHgc,
                         T &B, T &L, T &H, T &dB, T &dL, T &dH,
                         const TEarthDatum &earth=CK42);




//------------------------------------------------------------------------------
/*!
 * \brief Преобразование из Геодезической СК в ПГСК
 * \param B широта, [рад]
 * \param L долгота, [рад]
 * \param H высота, [м]
 * \param X координаты в ПГСК, [м]
 * \param Y координаты в ПГСК, [м]
 * \param Z координаты в ПГСК, [м]
 * \param earth параметры земли (по умолчанию СК-42)
 */
bool geodez2pgsk(const T &B,  const T &L,  const T &H,
                 T &X, T &Y, T &Z,
                 const TEarthDatum &earth=CK42);
/*!
 * \brief Преобразование из Геодезической СК в ПГСК
 * \param B широта, [рад]
 * \param L долгота, [рад]
 * \param H высота, [м]
 * \param dB скорость по широте, [рад/с]
 * \param dL скорость по долготе, [рад/с]
 * \param dH скорость по высоте, [м/с]
 * \param X координаты в ПГСК, [м]
 * \param Y координаты в ПГСК, [м]
 * \param Z координаты в ПГСК, [м]
 * \param VX скорости в ПГСК, [м\с]
 * \param VY скорости в ПГСК, [м\с]
 * \param VZ скорости в ПГСК, [м\с]
 * \param earth Параметры Земли
 */
bool geodez2pgsk(const T &B,  const T &L,  const T &H,
                 const T &dB, const T &dL, const T &dH,
                 T &X, T &Y, T &Z, T &VX, T &VY, T &VZ,
                 const TEarthDatum &earth=CK42);
//------------------------------------------------------------------------------
/*!
 * \brief Преобразование из ПГСК в Геодезическую СК
 * \param X координаты в ПГСК, [м]
 * \param Y координаты в ПГСК, [м]
 * \param Z координаты в ПГСК, [м]
 * \param B широта, [рад]
 * \param L долгота, [рад]
 * \param H высота, [м]
 * \param earth параметры земли (по умолчанию СК-42)
 */
bool pgsk2geodez(const T &X, const T &Y, const T &Z,
                 T &B, T &L, T &H, const TEarthDatum &earth=CK42);
/*!
 * \brief Преобразование из ПГСК в Геодезическую СК
 * \param X координаты в ПГСК, [м]
 * \param Y координаты в ПГСК, [м]
 * \param Z координаты в ПГСК, [м]
 * \param VX скорости в ПГСК, [м\с]
 * \param VY скорости в ПГСК, [м\с]
 * \param VZ скорости в ПГСК, [м\с]
 * \param B широта, [рад]
 * \param L долгота, [рад]
 * \param H высота, [м]
 * \param dB скорость по широте, [рад/с]
 * \param dL скорость по долготе, [рад/с]
 * \param dH скорость по высоте, [м/с]
 * \param earth параметры земли (по умолчанию СК-42)
 */
bool pgsk2geodez(const T &X, const T &Y, const T &Z,
                 const T &VX, const T &VY, const T &VZ,
                 T &B, T &L, T &H, T &dB, T &dL, T &dH,
                 const TEarthDatum &earth=CK42);




//------------------------------------------------------------------------------
/*!
 * \brief Преобразование из РТСК в МПСК
 * \param R дальность, [м]
 * \param B азимут, [рад]
 * \param E угол места, [рад]
 * \param dR производные сферических РТСК координат
 * \param dB производные сферических РТСК координат
 * \param dE производные сферических РТСК координат
 * \param X прямоуг. МПСК координаты, [м]
 * \param Y прямоуг. МПСК координаты, [м]
 * \param Z прямоуг. МПСК координаты, [м]
 * \param VX производные прямоуг. МПСК координат, [м/с]
 * \param VY производные прямоуг. МПСК координат, [м/с]
 * \param VZ производные прямоуг. МПСК координат, [м/с]
 */
bool rtsk2mpsk(const T &R,  const T &B,  const T &E,
               const T &dR, const T &dB, const T &dE,
               T &X, T &Y, T &Z, T &VX, T &VY, T &VZ);
//------------------------------------------------------------------------------                              |
/*!
 * \brief Преобразование из МПСК в РТСК
 * \param X прямоуг. МПСК координаты, [м]
 * \param Y прямоуг. МПСК координаты, [м]
 * \param Z прямоуг. МПСК координаты, [м]
 * \param VX производные прямоуг. МПСК координат, [м/с]
 * \param VY производные прямоуг. МПСК координат, [м/с]
 * \param VZ производные прямоуг. МПСК координат, [м/с]
 * \param R дальность, [м]
 * \param B азимут, [рад]
 * \param E угол места, [рад]
 * \param dR производные сферических РТСК координат
 * \param dB производные сферических РТСК координат
 * \param dE производные сферических РТСК координат
 */
bool mpsk2rtsk(const T &X, const T &Y, const T &Z,
               const T &VX, const T &VY, const T &VZ,
               T &R,  T &B,  T &E, T &dR, T &dB, T &dE);



//------------------------------------------------------------------------------
/*!
 * \brief расчёт матрицы поворота для перехода от ПГСК к МПСК
 * \param B геодезическая широта, [рад]
 * \param L геодезическая долгота, [рад]
 * \param Q азимут РСН, [рад]
 * \param e11 элементы матрицы поворота
 * \param e21 элементы матрицы поворота
 * \param e31 элементы матрицы поворота
 * \param e12 элементы матрицы поворота
 * \param e22 элементы матрицы поворота
 * \param e32 элементы матрицы поворота
 * \param e13 элементы матрицы поворота
 * \param e23 элементы матрицы поворота
 * \param e33 элементы матрицы поворота
 */
bool pgsk2mpskRotationMatrix(const T &B, const T &L, const T &Q,
                             T &e11, T &e21, T &e31,
                             T &e12, T &e22, T &e32,
                             T &e13, T &e23, T &e33);

//------------------------------------------------------------------------------
/*!
 * \brief расчёт матрицы перехода от ПГСК к МПСК
 * \param B геодезическая широта, [рад]
 * \param L геодезическая долгота, [рад]
 * \param H высота центра МПСК, [м]
 * \param Q азимут РСН, [рад]
 * \param e11 элементы матрицы перехода
 * \param e21 элементы матрицы перехода
 * \param e31 элементы матрицы перехода
 * \param e12 элементы матрицы перехода
 * \param e22 элементы матрицы перехода
 * \param e32 элементы матрицы перехода
 * \param e13 элементы матрицы перехода
 * \param e23 элементы матрицы перехода
 * \param e33 элементы матрицы перехода
 * \param shX смещение центра координат МПСК
 * \param shY смещение центра координат МПСК
 * \param shZ смещение центра координат МПСК
 * \param earth параметры земли
 */
bool pgsk2mpskMatr(const T &B, const T &L, const T &H, const T &Q,
                   T &e11, T &e21, T &e31,
                   T &e12, T &e22, T &e32,
                   T &e13, T &e23, T &e33,
                   T &shX, T &shY, T &shZ,
                   const TEarthDatum &earth=CK42);


//------------------------------------------------------------------------------
/*!
 * \brief Преобразование из ПГСК в МПСК
 * \param iX ПГСК координаты [м]
 * \param iY ПГСК координаты [м]
 * \param iZ ПГСК координаты [м]
 * \param shX смещение центра координат МПСК
 * \param shY смещение центра координат МПСК
 * \param shZ смещение центра координат МПСК
 * \param e11 элементы матрицы перехода
 * \param e21 элементы матрицы перехода
 * \param e31 элементы матрицы перехода
 * \param e12 элементы матрицы перехода
 * \param e22 элементы матрицы перехода
 * \param e32 элементы матрицы перехода
 * \param e13 элементы матрицы перехода
 * \param e23 элементы матрицы перехода
 * \param e33 элементы матрицы перехода
 * \param oX МПСК координаты [м]
 * \param oY МПСК координаты [м]
 * \param oZ МПСК координаты [м]
 */
bool pgsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ);
/*!
 * \brief Преобразование из ПГСК в МПСК
 * \param iX ПГСК координаты [м]
 * \param iY ПГСК координаты [м]
 * \param iZ ПГСК координаты [м]
 * \param iVX производные ПГСК координат [м/с]
 * \param iVY производные ПГСК координат [м/с]
 * \param iVZ производные ПГСК координат [м/с]
 * \param shX смещение центра координат МПСК
 * \param shY смещение центра координат МПСК
 * \param shZ смещение центра координат МПСК
 * \param e11 элементы матрицы перехода
 * \param e21 элементы матрицы перехода
 * \param e31 элементы матрицы перехода
 * \param e12 элементы матрицы перехода
 * \param e22 элементы матрицы перехода
 * \param e32 элементы матрицы перехода
 * \param e13 элементы матрицы перехода
 * \param e23 элементы матрицы перехода
 * \param e33 элементы матрицы перехода
 * \param oX МПСК координаты [м]
 * \param oY МПСК координаты [м]
 * \param oZ МПСК координаты [м]
 * \param oVX производные МПСК координат [м/с]
 * \param oVY производные МПСК координат [м/с]
 * \param oVZ производные МПСК координат [м/с]
 */
bool pgsk2mpsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ);

//------------------------------------------------------------------------------
/*!
 * \brief Преобразование из МПСК в ПГСК
 * \param iX ПГСК координаты [м]
 * \param iY ПГСК координаты [м]
 * \param iZ ПГСК координаты [м]
 * \param shX смещение центра координат МПСК
 * \param shY смещение центра координат МПСК
 * \param shZ смещение центра координат МПСК
 * \param e11 элементы матрицы перехода
 * \param e21 элементы матрицы перехода
 * \param e31 элементы матрицы перехода
 * \param e12 элементы матрицы перехода
 * \param e22 элементы матрицы перехода
 * \param e32 элементы матрицы перехода
 * \param e13 элементы матрицы перехода
 * \param e23 элементы матрицы перехода
 * \param e33 элементы матрицы перехода
 * \param oX МПСК координаты [м]
 * \param oY МПСК координаты [м]
 * \param oZ МПСК координаты [м]
 */
bool mpsk2pgsk(const T &iX,  const T &iY,  const T &iZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ);
/*!
 * \brief Преобразование из МПСК в ПГСК
 * \param iX ПГСК координаты [м]
 * \param iY ПГСК координаты [м]
 * \param iZ ПГСК координаты [м]
 * \param iVX производные ПГСК координат [м/с]
 * \param iVY производные ПГСК координат [м/с]
 * \param iVZ производные ПГСК координат [м/с]
 * \param shX смещение центра координат МПСК
 * \param shY смещение центра координат МПСК
 * \param shZ смещение центра координат МПСК
 * \param e11 элементы матрицы перехода
 * \param e21 элементы матрицы перехода
 * \param e31 элементы матрицы перехода
 * \param e12 элементы матрицы перехода
 * \param e22 элементы матрицы перехода
 * \param e32 элементы матрицы перехода
 * \param e13 элементы матрицы перехода
 * \param e23 элементы матрицы перехода
 * \param e33 элементы матрицы перехода
 * \param oX МПСК координаты [м]
 * \param oY МПСК координаты [м]
 * \param oZ МПСК координаты [м]
 * \param oVX производные МПСК координат [м/с]
 * \param oVY производные МПСК координат [м/с]
 * \param oVZ производные МПСК координат [м/с]
 */
bool mpsk2pgsk(const T &iX,  const T &iY,  const T &iZ,
               const T &iVX, const T &iVY, const T &iVZ,
               const T &shX, const T &shY, const T &shZ,
               const T &e11, const T &e21, const T &e31,
               const T &e12, const T &e22, const T &e32,
               const T &e13, const T &e23, const T &e33,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ);


//------------------------------------------------------------------------------
/*!
 * \brief Определение дальности по Земле от точки X до точки Y
 * \param Bx геодезические координаты точки Х
 * \param Lx геодезические координаты точки Х
 * \param Hx геодезические координаты точки Х
 * \param By геодезические координаты точки Y
 * \param Ly геодезические координаты точки Y
 * \param Hy геодезические координаты точки Y
 * \param earth параметры земли
 */
T sphereDistance(const T &Bx, const T &Lx, const T &Hx,
                 const T &By, const T &Ly, const T &Hy,
                 const TEarthDatum &earth=CK42);


/*!
 * \brief Расчёт матрицы перехода от МПСК <=> СТСК
 * \param mB Геодезическая широта начала МПСК
 * \param mL Геодезическая долгота начала МПСК
 * \param mH Геодезическая высота начала МПСК
 * \param sB Геодезическая широта начала СТСК
 * \param sL Геодезическая долгота начала СТСК
 * \param sH Геодезическая высота начала СТСК
 * \param sQ Азимут стрельбы
 * \param e11 элемент матрицы поворота
 * \param e21 элемент матрицы поворота
 * \param e31 элемент матрицы поворота
 * \param e12 элемент матрицы поворота
 * \param e22 элемент матрицы поворота
 * \param e32 элемент матрицы поворота
 * \param e13 элемент матрицы поворота
 * \param e23 элемент матрицы поворота
 * \param e33 элемент матрицы поворота
 * \param shX элемент вектора смещения
 * \param shY элемент вектора смещения
 * \param shZ элемент вектора смещения
 * \param earth Параметры Земли
 */
bool mpsk2stskMatr(const T &mB, const T &mL, const T &mH,
                   const T &sB, const T &sL, const T &sH, const T &sQ,
                   T &e11, T &e21, T &e31,
                   T &e12, T &e22, T &e32,
                   T &e13, T &e23, T &e33,
                   T &shX, T &shY, T &shZ, const TEarthDatum &earth=CK42);

/*!
 * \brief Преобразование МПСК <=> СТСК (альтернативное преобразование)
 * \param iX координата X в МПСК
 * \param iY координата Y в МПСК
 * \param iZ координата Z в МПСК
 * \param shX смещение центра координат СТСК по X
 * \param shY смещение центра координат СТСК по Y
 * \param shZ смещение центра координат СТСК по Z
 * \param e11 элемент матрицы поворота
 * \param e21 элемент матрицы поворота
 * \param e31 элемент матрицы поворота
 * \param e12 элемент матрицы поворота
 * \param e22 элемент матрицы поворота
 * \param e32 элемент матрицы поворота
 * \param e13 элемент матрицы поворота
 * \param e23 элемент матрицы поворота
 * \param e33 элемент матрицы поворота
 * \param oX координата X в СТСК
 * \param oY координата Y в СТСК
 * \param oZ координата Z в СТСК
 * \return Флаг результата
 */
bool mpsk2stskA(const T &iX,  const T &iY,  const T &iZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ);
/*!
 * \brief Преобразование МПСК <=> СТСК (альтернативное преобразование)
 * \param iX координата X в МПСК
 * \param iY координата Y в МПСК
 * \param iZ координата Z в МПСК
 * \param iVX составляющая скорости по X в МПСК
 * \param iVY составляющая скорости по Y в МПСК
 * \param iVZ составляющая скорости по Z в МПСК
 * \param shX смещение центра координат СТСК по X
 * \param shY смещение центра координат СТСК по Y
 * \param shZ смещение центра координат СТСК по Z
 * \param e11 элемент матрицы поворота
 * \param e21 элемент матрицы поворота
 * \param e31 элемент матрицы поворота
 * \param e12 элемент матрицы поворота
 * \param e22 элемент матрицы поворота
 * \param e32 элемент матрицы поворота
 * \param e13 элемент матрицы поворота
 * \param e23 элемент матрицы поворота
 * \param e33 элемент матрицы поворота
 * \param oX координата X в СТСК
 * \param oY координата Y в СТСК
 * \param oZ координата Z в СТСК
 * \param oVX составляющая скорости по X в СТСК
 * \param oVY составляющая скорости по Y в СТСК
 * \param oVZ составляющая скорости по Z в СТСК
 * \return Флаг результата
 */
bool mpsk2stskA(const T &iX,  const T &iY,  const T &iZ,
                const T &iVX, const T &iVY, const T &iVZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ);
/*!
 * \brief Преобразование МПСК <=> СТСК (альтернативное преобразование)
 * \param iX координата X в МПСК
 * \param iY координата Y в МПСК
 * \param iZ координата Z в МПСК
 * \param shX смещение центра координат СТСК по X
 * \param shY смещение центра координат СТСК по Y
 * \param shZ смещение центра координат СТСК по Z
 * \param e11 элемент матрицы поворота
 * \param e21 элемент матрицы поворота
 * \param e31 элемент матрицы поворота
 * \param e12 элемент матрицы поворота
 * \param e22 элемент матрицы поворота
 * \param e32 элемент матрицы поворота
 * \param e13 элемент матрицы поворота
 * \param e23 элемент матрицы поворота
 * \param e33 элемент матрицы поворота
 * \param oX координата X в СТСК
 * \param oY координата Y в СТСК
 * \param oZ координата Z в СТСК
 * \return Флаг результата
 */
bool stsk2mpskA(const T &iX,  const T &iY,  const T &iZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ);
/*!
 * \brief Преобразование МПСК <=> СТСК (альтернативное преобразование)
 * \param iX координата X в МПСК
 * \param iY координата Y в МПСК
 * \param iZ координата Z в МПСК
 * \param iVX составляющая скорости по X в МПСК
 * \param iVY составляющая скорости по Y в МПСК
 * \param iVZ составляющая скорости по Z в МПСК
 * \param shX смещение центра координат СТСК по X
 * \param shY смещение центра координат СТСК по Y
 * \param shZ смещение центра координат СТСК по Z
 * \param e11 элемент матрицы поворота
 * \param e21 элемент матрицы поворота
 * \param e31 элемент матрицы поворота
 * \param e12 элемент матрицы поворота
 * \param e22 элемент матрицы поворота
 * \param e32 элемент матрицы поворота
 * \param e13 элемент матрицы поворота
 * \param e23 элемент матрицы поворота
 * \param e33 элемент матрицы поворота
 * \param oX координата X в СТСК
 * \param oY координата Y в СТСК
 * \param oZ координата Z в СТСК
 * \param oVX составляющая скорости по X в СТСК
 * \param oVY составляющая скорости по Y в СТСК
 * \param oVZ составляющая скорости по Z в СТСК
 * \return Флаг результата
 */
bool stsk2mpskA(const T &iX,  const T &iY,  const T &iZ,
                const T &iVX, const T &iVY, const T &iVZ,
                const T &shX, const T &shY, const T &shZ,
                const T &e11, const T &e21, const T &e31,
                const T &e12, const T &e22, const T &e32,
                const T &e13, const T &e23, const T &e33,
                T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ);



//------------------------------------------------------------------------------
/*!
 * \brief Преобразование координат из РТСК в АПСК
 * \param r [м]
 * \param b [рад]
 * \param e [рад]
 * \param dR [м\с]
 * \param dB [рад\с]
 * \param dE [рад\с]
 * \param x [м]
 * \param y [м]
 * \param z [м]
 * \param VX [м\с]
 * \param VY [м\с]
 * \param VZ [м\с]
 */
bool rtsk2apsk(const T &r, const T &b, const T &e,
               const T &dR, const T &dB, const T &dE,
               T &x, T &y, T &z, T &VX, T &VY, T &VZ);

//------------------------------------------------------------------------------
/*!
 * \brief Преобразование координат из АПСК в РТСК
 * \param x [м]
 * \param y [м]
 * \param z [м]
 * \param VX [м\с]
 * \param VY [м\с]
 * \param VZ [м\с]
 * \param r [м]
 * \param b [рад]
 * \param e [рад]
 * \param dR [м\с]
 * \param dB [рад\с]
 * \param dE [рад\с]
 */
bool apsk2rtsk(const T &x, const T &y,  const T &z,
               const T &VX, const T &VY, const T &VZ,
               T &r, T &b, T &e, T &dR, T &dB, T &dE);


//------------------------------------------------------------------------------
/*!
 * \brief Преобразование координат из АПСК в ОБСК
 * \param x [м]
 * \param y [м]
 * \param z [м]
 * \param VX [м\с]
 * \param VY [м\с]
 * \param VZ [м\с]
 * \param oR координаты ОБСК [м]
 * \param oU координаты ОБСК [м]
 * \param oV координаты ОБСК [м]
 * \param odR производные координат ОБСК [м/c]
 * \param odU производные координат ОБСК [м/c]
 * \param odV производные координат ОБСК [м/c]
 */
bool apsk2obsk(const T &x, const T &y,  const T &z,
               const T &VX, const T &VY, const T &VZ,
               T &oR, T &oU, T &oV, T &odR, T &odU, T &odV);


//------------------------------------------------------------------------------
/*!
 * \brief Преобразование координат из АПСК в ОБСК
 * \param iR координаты ОБСК [м]
 * \param iU координаты ОБСК [м]
 * \param iV координаты ОБСК [м]
 * \param idR производные координат ОБСК [м/c]
 * \param idU производные координат ОБСК [м/c]
 * \param idV производные координат ОБСК [м/c]
 * \param oX [м]
 * \param oY [м]
 * \param oZ [м]
 * \param oVX [м\с]
 * \param oVY [м\с]
 * \param oVZ [м\с]
 */
bool obsk2apsk(const T &iR, const T &iU,  const T &iV,
               const T &idR, const T &idU, const T &idV,
               T &oX, T &oY, T &oZ, T &oVX, T &oVY, T &oVZ);




//------------------------------------------------------------------------------
/*!
 * \brief Расчёт поправочных значений координат с учётом разницы эпох
 * \param imonth Месяц текущей эпохи
 * \param iyear Год текущей эпохи
 * \param omonth Месяц искомой эпохи
 * \param oyear Год искомой эпохи
 * \param dX Координата Х в искомой эпохе
 * \param dY Координата Y в искомой эпохе
 * \param dZ Координата Z в искомой эпохе
 * \param earth Параметры Земли
 * \return Применялась ли операция коррекции (или координаты прежние)
 */
bool calcPGSKDeviationsOnEpoch(ruint8 imonth, ruint16 iyear,
                               ruint8 omonth, ruint16 oyear,
                               T &dX, T &dY,  T &dZ, const TEarthDatum &earth);



//------------------------------------------------------------------------------
/*!
 * \brief Формирование матриц коррекции координат для переходов между датумами
 * \param fromEarth Параметры Земли для исходных координат
 * \param toEarth Параметры Земли для искомых координат
 * \param e11 Элементы матрицы поворота
 * \param e21 Элементы матрицы поворота
 * \param e31 Элементы матрицы поворота
 * \param e12 Элементы матрицы поворота
 * \param e22 Элементы матрицы поворота
 * \param e32 Элементы матрицы поворота
 * \param e13 Элементы матрицы поворота
 * \param e23 Элементы матрицы поворота
 * \param e33 Элементы матрицы поворота
 * \param shX Элементы матрицы сдвига
 * \param shY Элементы матрицы сдвига
 * \param shZ Элементы матрицы сдвига
 * \return Результат выполнения операции: 0 - неудача; >0 - ок;
 *                                       <0 - операция не имеет смысла
 */
rint32 earthDevMatrix(const TEarthDatum &fromEarth, const TEarthDatum &toEarth,
                      T &e11, T &e21, T &e31,
                      T &e12, T &e22, T &e32,
                      T &e13, T &e23, T &e33,
                      T &shX, T &shY, T &shZ);
/*!
 * \brief Расчёт координат с поправками датумов
 * \param iX Координата Х в текущем датуме
 * \param iY Координата Y в текущем датуме
 * \param iZ Координата Z в текущем датуме
 * \param oX Координата Х в искомом датуме
 * \param oY Координата Y в искомом датуме
 * \param oZ Координата Z в искомом датуме
 * \param e11 Элементы матрицы поворота
 * \param e21 Элементы матрицы поворота
 * \param e31 Элементы матрицы поворота
 * \param e12 Элементы матрицы поворота
 * \param e22 Элементы матрицы поворота
 * \param e32 Элементы матрицы поворота
 * \param e13 Элементы матрицы поворота
 * \param e23 Элементы матрицы поворота
 * \param e33 Элементы матрицы поворота
 * \param shX Элементы матрицы сдвига
 * \param shY Элементы матрицы сдвига
 * \param shZ Элементы матрицы сдвига
 */
void calculateDeviations(const T &iX, const T &iY, const T &iZ,
                         T &oX, T &oY, T &oZ,
                         T &e11, T &e21, T &e31,
                         T &e12, T &e22, T &e32,
                         T &e13, T &e23, T &e33,
                         T &shX, T &shY, T &shZ);



} //namespace RGeoMath

} //namespace RMath


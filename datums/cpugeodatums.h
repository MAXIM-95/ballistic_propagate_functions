#pragma once

#include "crdefs.h"
#include "math/cpumath.h"

using namespace RMath;

namespace RModels {

namespace RDatums {

/*!
 * \brief Параметры Земли
 */
class TEarthDatum {
public:
  /*!
   * \brief Типы параметров Земли
   */
  enum CEarthTypes {
    ET_CK42    = 0, ///<СК-42
    ET_PZ90    = 1, ///<ПЗ-90
    ET_PZ90_02 = 2, ///<ПЗ-90.02
    ET_PZ90_11 = 3, ///<ПЗ-90.11
    ET_GSK2011 = 4, ///<ГСК-2011
    ET_WGS84   = 5, ///<WGS-84
    ET_USERTYPES    ///<Пользовательские
  };
  
  /*!
   * \brief Трансформация координат между параметрами Земли
   * \param ix Координата X положения для старых ПЗ
   * \param iy Координата Y положения для старых ПЗ
   * \param iz Координата Z положения для старых ПЗ
   * \param ifrom Описание старых ПЗ
   * \param ito Описание новых ПЗ
   * \return Признак результата операции
   *         (результат выводится во входные параметры)
   */
  static bool transformByEarthDatums(T &ix, T &iy, T &iz,
                                     const TEarthDatum::CEarthTypes ifrom,
                                     const TEarthDatum::CEarthTypes ito);
  TEarthDatum(rint32 itp=ET_CK42);
  /*!
   * \brief Создание параметров Земли
   * \param itp Тип параметров
   * \param ij2 Коэффициент 2-й зональной гармоники
   * \param iea Коэффициент сжатия Земли
   * \param iggc Геоцентрическая гравитационная постоянная Земли
   * \param ibsa Большая полуось(экваториальный радиус)
   * \param iave Угловая скорость вращения Земли
   */
  TEarthDatum(rint32 itp, const T &ij2, const T &iea, const T &iggc, const T &ibsa,
              const T &iave);
  TEarthDatum(const TEarthDatum& copy);
  TEarthDatum & operator=(const TEarthDatum &copy);

  inline friend bool operator==(const TEarthDatum &a, const TEarthDatum &b)
  {
    if ( (a._type<ET_USERTYPES) && (a._type==b._type) ) return true;

    return ( (a._GeocentricalGravityConst == b._GeocentricalGravityConst) &&
             (a._EllipticityEarth         == b._EllipticityEarth) &&
             (a._BigSemiAxis              == b._BigSemiAxis) &&
             (a._AngularVelocityEarth     == b._AngularVelocityEarth) &&
             (a._J2                       == b._J2) );
  }
  inline friend bool operator!=(const TEarthDatum &a, const TEarthDatum &b)
  {
    return !operator ==(a, b);
  }

  /*!
   * \brief Тип параметров Земли
   * \return Тип параметров Земли
   */
  rint32 type()const;

  /*!
   * \brief Коэффициент сжатия Земли
   * \return Коэффициент сжатия Земли
   */
  T EllipticityEarth()const;

  /*!
   * \brief Производная коэффициента сжатия Земли
   * \return Производная коэффициента сжатия Земли
   */ 
  T EllipticityEarthDeriv()const;

  /*!
   * \brief Коэффициент 2-й зональной гармоники
   * \return Коэффициент 2-й зональной гармоники
   */
  T J2()const;

  /*!
   * \brief Масса Земли(включая массу атмосферы)
   * \return Масса Земли(включая массу атмосферы)
   */
  T EarthWeight()const;

  /*!
   * \brief Большая полуось(экваториальный радиус)
   * \return Большая полуось(экваториальный радиус)
   */
  T BigSemiAxis()const;

  /*!
   * \brief Малая полуось(полярный радиус)
   * \return Малая полуось(полярный радиус)
   */
  T SmallSemiAxis()const;

  /*!
   * \brief Коэффициент несферичности Земли
   * \return Коэффициент несферичности Земли
   */
  T Nonsphericity()const;

  /*!
   * \brief Эксцентриситет
   * \return Эксцентриситет
   */
  T Eccentricity()const;

  /*!
   * \brief Геоцентрическая гравитационная постоянная Земли
   * \return Геоцентрическая гравитационная постоянная Земли
   */
  T GravityConstant()const;

  /*!
   * \brief Параметр для учёта влияния гравитационного поля Земли(а00)
   * \return Параметр для учёта влияния гравитационного поля Земли(а00)
   */
  T a00()const;
  /*!
   * \brief Параметр для учёта влияния гравитационного поля Земли(а20)
   * \return Параметр для учёта влияния гравитационного поля Земли(а20)
   */
  T a20()const;


  /*!
   * \brief Угловая скорость вращения Земли
   * \return Угловая скорость вращения Земли
   */
  T AngularVelocity()const;

  /*!
   * \brief Квадрат угловой скорости вращения Земли
   * \return Квадрат угловой скорости вращения Земли
   */
  T AngularVelocity2()const;

  /*!
   * \brief Средний радиус Земли
   * \return Средний радиус Земли
   */
  static T meanRadius();

  /*!
   * \brief Отношение массы Земли к массе Луны
   * \return Отношение массы Земли к массе Луны
   */
  static T MoonMassKoef();

  /*!
   * \brief Ускорение свободного падения на Земле
   * \return Ускорение свободного падения на Земле
   */
  static T G0();

  /*!
   * \brief Проверка корректности создания объекта
   * \return Признак корректности создания
   */
  bool isValid()const;

private:
  rint32 _type;
  bool   _valid;
  T      _a00,                       ///<параметры для учёта гравитации
         _a20,                       ///<параметры для учёта гравитации
         _GeocentricalGravityConst,  ///<геоцентрическая гравитационная постоянная Земли
         //масса Земли(вкл. массу атмосферы) умн. на грав. пост.[м^3\с^2]
         _EllipticityEarth,          ///<коэффициент сжатия Земли
         _BigSemiAxis,               ///<большая полуось(экваториальный радиус),[м]
         _AngularVelocityEarth,      ///<угловая скорость вращения Земли, [рад/с]
         _J2,                        ///<коэффициент 2-й зональной гармоники
         _EarthWeight,               ///<масса Земли(включая массу атмосферы)
         _EllipticityEarthDeriv,     ///<производная коэффициента сжатия Земли
         _SmallSemiAxis,             ///<малая полуось, [м]
         _Eccentricity,              ///<эксцентриситет
         _Nonsphericity,             ///<Коэффициент несферичности Земли
         _AngularVelocityEarth2;     ///<квадрат угловой скорости вращения Земли, [рад^2/с^2]

};

  extern TEarthDatum CK42;
  extern TEarthDatum PZ90;
  extern TEarthDatum PZ90_02;
  extern TEarthDatum PZ90_11;
  extern TEarthDatum GSK2011;
  extern TEarthDatum WGS84;


} //namespace RDatums

} //namespace RModels


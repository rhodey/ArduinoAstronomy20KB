#ifndef __ASTRONOMY_ENGINE_H__
#define __ASTRONOMY_ENGINE_H__

#if (ARDUINO >= 100)
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

#include <SPIFlash.h>

/**
 * @brief Indicates success/failure of an Astronomy Engine function call.
 */
typedef enum
{
    ASTRO_SUCCESS,                  /**< The operation was successful. */
    ASTRO_NOT_INITIALIZED,          /**< A placeholder that can be used for data that is not yet initialized. */
    ASTRO_INVALID_BODY,             /**< The celestial body was not valid. Different sets of bodies are supported depending on the function. */
    ASTRO_NO_CONVERGE,              /**< A numeric solver failed to converge. This should not happen unless there is a bug in Astronomy Engine. */
    ASTRO_BAD_TIME,                 /**< The provided date/time is outside the range allowed by this function. */
    ASTRO_BAD_VECTOR,               /**< Vector magnitude is too small to be normalized into a unit vector. */
    ASTRO_SEARCH_FAILURE,           /**< Search was not able to find an ascending root crossing of the function in the specified time interval. */
    ASTRO_EARTH_NOT_ALLOWED,        /**< The Earth cannot be treated as a celestial body seen from an observer on the Earth itself. */
    ASTRO_NO_MOON_QUARTER,          /**< No lunar quarter occurs inside the specified time range. */
    ASTRO_WRONG_MOON_QUARTER,       /**< Internal error: Astronomy_NextMoonQuarter found the wrong moon quarter. */
    ASTRO_INTERNAL_ERROR,           /**< A self-check failed inside the code somewhere, indicating a bug needs to be fixed. */
    ASTRO_INVALID_PARAMETER,        /**< A parameter value passed to a function was not valid. */
    ASTRO_FAIL_APSIS,               /**< Special-case logic for finding Neptune/Pluto apsis failed. */
    ASTRO_BUFFER_TOO_SMALL,         /**< A provided buffer's size is too small to receive the requested data. */
    ASTRO_OUT_OF_MEMORY             /**< An attempt to allocate memory failed. */
}
astro_status_t;

/**
 * @brief A date and time used for astronomical calculations.
 *
 * This type is of fundamental importance to Astronomy Engine.
 * It is used to represent dates and times for all astronomical calculations.
 * It is also included in the values returned by many Astronomy Engine functions.
 *
 * To create a valid astro_time_t value from scratch, call #Astronomy_MakeTime
 * (for a given calendar date and time) or #Astronomy_CurrentTime (for the system's
 * current date and time).
 *
 * To adjust an existing astro_time_t by a certain real number of days,
 * call #Astronomy_AddDays.
 *
 * The astro_time_t type contains `ut` to represent Universal Time (UT1/UTC) and
 * `tt` to represent Terrestrial Time (TT, also known as *ephemeris time*).
 * The difference `tt-ut` is known as *&Delta;T*, and is obtained from
 * a model provided by the
 * [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).
 *
 * Both `tt` and `ut` are necessary for performing different astronomical calculations.
 * Indeed, certain calculations (such as rise/set times) require both time scales.
 * See the documentation for the `ut` and `tt` fields for more detailed information.
 *
 * In cases where astro_time_t is included in a structure returned by
 * a function that can fail, the astro_status_t field `status` will contain a value
 * other than `ASTRO_SUCCESS`; in that case the `ut` and `tt` will hold `NAN` (not a number).
 * In general, when there is an error code stored in a struct field `status`, the
 * caller should ignore all other values in that structure, including the `ut` and `tt`
 * inside astro_time_t.
 */
typedef struct
{
    /**
     * @brief   UT1/UTC number of days since noon on January 1, 2000.
     *
     * The floating point number of days of Universal Time since noon UTC January 1, 2000.
     * Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
     * not exactly equivalent; UTC and UT1 can disagree by up to &plusmn;0.9 seconds.
     * This approximation is sufficient for the accuracy requirements of Astronomy Engine.
     *
     * Universal Time Coordinate (UTC) is the international standard for legal and civil
     * timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
     * UTC is kept in sync with unpredictable observed changes in the Earth's rotation
     * by occasionally adding leap seconds as needed.
     *
     * UT1 is an idealized time scale based on observed rotation of the Earth, which
     * gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
     * large scale weather events like hurricanes, and internal seismic and convection effects.
     * Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
     * is adjusted by a scheduled whole number of leap seconds as needed.
     *
     * The value in `ut` is appropriate for any calculation involving the Earth's rotation,
     * such as calculating rise/set times, culumination, and anything involving apparent
     * sidereal time.
     *
     * Before the era of atomic timekeeping, days based on the Earth's rotation
     * were often known as *mean solar days*.
     */
    double ut;

    /**
     * @brief   Terrestrial Time days since noon on January 1, 2000.
     *
     * Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
     * In this system, days are not based on Earth rotations, but instead by
     * the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
     * divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
     * for changes in the Earth's rotation.
     *
     * The value in `tt` is used for calculations of movements not involving the Earth's rotation,
     * such as the orbits of planets around the Sun, or the Moon around the Earth.
     *
     * Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
     */
    double tt;

    /**
     * @brief   For internal use only. Used to optimize Earth tilt calculations.
     */
    double psi;

    /**
     * @brief   For internal use only.  Used to optimize Earth tilt calculations.
     */
    double eps;
}
astro_time_t;

/**
 * @brief A celestial body.
 */
typedef enum
{
    BODY_INVALID = -1,      /**< An invalid or undefined celestial body. */
    BODY_MERCURY,           /**< Mercury */
    BODY_VENUS,             /**< Venus */
    BODY_EARTH,             /**< Earth */
    BODY_MARS,              /**< Mars */
    BODY_JUPITER,           /**< Jupiter */
    BODY_SATURN,            /**< Saturn */
    BODY_URANUS,            /**< Uranus */
    BODY_NEPTUNE,           /**< Neptune */
    BODY_PLUTO,             /**< Pluto */
    BODY_SUN,               /**< Sun */
    BODY_MOON,              /**< Moon */
    BODY_EMB,               /**< Earth/Moon Barycenter */
    BODY_SSB                /**< Solar System Barycenter */
}
astro_body_t;

/**
 * @brief The location of an observer on (or near) the surface of the Earth.
 *
 * This structure is passed to functions that calculate phenomena as observed
 * from a particular place on the Earth.
 *
 * You can create this structure directly, or you can call the convenience function
 * #Astronomy_MakeObserver# to create one for you.
 */
typedef struct
{
    double latitude;        /**< Geographic latitude in degrees north (positive) or south (negative) of the equator. */
    double longitude;       /**< Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. */
    double height;          /**< The height above (positive) or below (negative) sea level, expressed in meters. */
}
astro_observer_t;

/**
 * @brief   Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.
 *
 * The Earth's equator is not always in the same plane due to precession and nutation.
 *
 * Sometimes it is useful to have a fixed plane of reference for equatorial coordinates
 * across different calendar dates.  In these cases, a fixed *epoch*, or reference time,
 * is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers
 * to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.
 *
 * For some other purposes, it is more helpful to represent coordinates using the Earth's
 * equator exactly as it is on that date. For example, when calculating rise/set times
 * or horizontal coordinates, it is most accurate to use the orientation of the Earth's
 * equator at that same date and time. For these uses, Astronomy Engine allows *of-date*
 * calculations.
 */
typedef enum
{
    EQUATOR_J2000,      /**< Represent equatorial coordinates in the J2000 epoch. */
    EQUATOR_OF_DATE     /**< Represent equatorial coordinates using the Earth's equator at the given date and time. */
}
astro_equator_date_t;

/**
 * @brief   Aberration calculation options.
 *
 * [Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect
 * causing the apparent direction of an observed body to be shifted due to transverse
 * movement of the Earth with respect to the rays of light coming from that body.
 * This angular correction can be anywhere from 0 to about 20 arcseconds,
 * depending on the position of the observed body relative to the instantaneous
 * velocity vector of the Earth.
 *
 * Some Astronomy Engine functions allow optional correction for aberration by
 * passing in a value of this enumerated type.
 *
 * Aberration correction is useful to improve accuracy of coordinates of
 * apparent locations of bodies seen from the Earth.
 * However, because aberration affects not only the observed body (such as a planet)
 * but the surrounding stars, aberration may be unhelpful (for example)
 * for determining exactly when a planet crosses from one constellation to another.
 */
typedef enum
{
    ABERRATION,     /**< Request correction for aberration. */
    NO_ABERRATION   /**< Do not correct for aberration. */
}
astro_aberration_t;

/**
 * @brief Equatorial angular coordinates.
 *
 * Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context),
 * oriented with respect to the projection of the Earth's equator onto the sky.
 */
typedef struct
{
    astro_status_t status;  /**< `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. */
    double ra;              /**< right ascension in sidereal hours. */
    double dec;             /**< declination in degrees */
    double dist;            /**< distance to the celestial body in AU. */
}
astro_equatorial_t;

/**
 * @brief A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
 */
typedef struct
{
    astro_status_t status;  /**< `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. */
    double x;               /**< The Cartesian x-coordinate of the vector in AU. */
    double y;               /**< The Cartesian y-coordinate of the vector in AU. */
    double z;               /**< The Cartesian z-coordinate of the vector in AU. */
    astro_time_t t;         /**< The date and time at which this vector is valid. */
}
astro_vector_t;

/**
 * @brief Contains a rotation matrix that can be used to transform one coordinate system to another.
 */
typedef struct
{
    astro_status_t status;  /**< `ASTRO_SUCCESS` if this struct is valid; otherwise an error code. */
    double rot[3][3];       /**< A normalized 3x3 rotation matrix. */
}
astro_rotation_t;

/**
 * @brief Selects whether to correct for atmospheric refraction, and if so, how.
 */
typedef enum
{
    REFRACTION_NONE,    /**< No atmospheric refraction correction (airless). */
    REFRACTION_NORMAL,  /**< Recommended correction for standard atmospheric refraction. */
    REFRACTION_JPLHOR   /**< Used only for compatibility testing with JPL Horizons online tool. */
}
astro_refraction_t;

/**
 * @brief Coordinates of a celestial body as seen by a topocentric observer.
 *
 * Contains horizontal and equatorial coordinates seen by an observer on or near
 * the surface of the Earth (a topocentric observer).
 * Optionally corrected for atmospheric refraction.
 */
typedef struct
{
    double azimuth;     /**< Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. */
    double altitude;    /**< Angle in degrees above (positive) or below (negative) the observer's horizon. */
    double ra;          /**< Right ascension in sidereal hours. */
    double dec;         /**< Declination in degrees. */
}
astro_horizon_t;

double Astronomy_VectorLength(astro_vector_t vector);
double Astronomy_Refraction(astro_refraction_t refraction, double altitude);
astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second);
astro_time_t Astronomy_AddDays(astro_time_t time, double days);
astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time);
astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time, astro_aberration_t aberration);
astro_vector_t Astronomy_GeoMoon(astro_time_t time);

astro_equatorial_t Astronomy_Equator(
    SPIFlash flash,
    astro_body_t body,
    astro_time_t *time,
    astro_observer_t observer,
    astro_equator_date_t equdate,
    astro_aberration_t aberration
);

astro_horizon_t Astronomy_Horizon(
    SPIFlash flash,
    astro_time_t *time,
    astro_observer_t observer,
    double ra,
    double dec,
    astro_refraction_t refraction
);

#define QQQ 6
void blinkLED(int pin, int delayMs);

#endif

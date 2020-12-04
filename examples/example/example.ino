#include <ArduinoAstronomy.h>
#include <SPIFlash.h>

#define BAUD 115200
#define SPI_CS 10
#define DEVICE_ID 0xEF40 // W25Q32FV

SPIFlash flash(SPI_CS, DEVICE_ID);

void setup() {
  Serial.begin(BAUD);
  if (!flash.initialize()) {
    Serial.print("unexpected device ID: ");
    Serial.println(flash.readDeviceId(), HEX);
  }
}

void horizon(astro_time_t time, astro_observer_t observer, astro_body_t body) {
  astro_equatorial_t equ_ofdate = Astronomy_Equator(flash, body, &time, observer, EQUATOR_OF_DATE, ABERRATION);
  astro_horizon_t horizon = Astronomy_Horizon(flash, &time, observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);

  Serial.print("(azimuth, altitude) => (");
  Serial.print(horizon.azimuth);
  Serial.print(", ");
  Serial.print(horizon.altitude);
  Serial.println(")");
}

/*
earth = 19130
sun   = 19670
vsop  = 20738 ... 64%
moon  = 27486 ... 85%
venus + merc = 30046 ... 93%
venus + sun  = 31040 ... 96%
venus + moon = 31280 ... 96%
*/
void loop(void) {
  astro_time_t time = Astronomy_MakeTime(2020, 12, 4, 11, 30, 0); // year, month, day, utc hours, minutes, seconds
  astro_observer_t observer;
  observer.height = 0.0;
  observer.latitude = 41.4857608;
  observer.longitude = -71.3102154;

  astro_body_t body[] = {
    BODY_SUN, BODY_MOON, BODY_MERCURY, BODY_VENUS, BODY_MARS,
    BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE/*, BODY_PLUTO*/
  };

  for (int i = 0; i < 9; i++) {
    horizon(time, observer, body[i]);
  }

  Serial.println("done.");
}

# ArduinoAstronomy20KB
This library is a fork of [Astronomy-Engine](https://github.com/cosinekitty/astronomy) which moves `10KB` of solar system model parameters onto SPI flash memory using LowPowerLab's [SPIFlash](https://github.com/LowPowerLab/SPIFlash) library. This is not a full port of [Astronomy-Engine C/C++](https://github.com/cosinekitty/astronomy/blob/master/source/c/README.md), presently [Pluto is missing](https://github.com/rhodey/ArduinoAstronomy/issues/1), however azimuth and altitude calculations for all 8 planets plus the Sun and Moon are supported. Program size including calculations for the seven planets plus the Sun and Moon is `31424B`, seven planets `28696B`, Moon only `27488B`, Sun only `19672B`.

## Download & Install
To install a library for the Arduino IDE you must know where the `libraries/` directory is located then copy or link into it. I have installed the Arduino IDE using snap on Ubuntu so:
```
$ git clone https://github.com/rhodey/ArduinoAstronomy
$ ln -s $(pwd)/ArduinoAstronomy ~/snap/arduino/current/Arduino/libraries/ArduinoAstronomy
```

## Prepare Flash
Open `examples/flash/flash.ino` in the Arduino IDE, modify `SPI_CS` and `DEVICE_ID` if necessary then connect your Arduino and upload the program. Next install [Node.js](https://nodejs.org/en/download/) and run the following commands correcting the Arduino serial port path if necessary:
```
$ npm install
$ node lib/js/flash.js /dev/ttyUSB0
```

## Example
The following example can also be found in `examples/example/example.ino`:
```c
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
```

## License
MIT

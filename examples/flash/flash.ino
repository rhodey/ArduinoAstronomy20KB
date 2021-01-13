#include <SPIFlash.h>

#define BAUD 115200
#define P_SPI_CS 1
#define DEVICE_ID 0xEF40

SPIFlash flash(P_SPI_CS, DEVICE_ID);

void setup() {
  pinMode(P_SPI_CS, OUTPUT);
  Serial.begin(BAUD);
  while (!Serial) { delay(10); }
  if (!flash.initialize()) {
    Serial.print("unexpected device ID: ");
    Serial.println(flash.readDeviceId(), HEX);
  }
}

bool writeInt(long address, int val) {
  union {
    int val;
    uint8_t byte[2];
  } read;
  const uint8_t* p = (const uint8_t*)(const void*)&val;
  flash.writeBytes(address, p, 2);
  flash.readBytes(address, read.byte, 2);
  return read.val == val;
}

bool writeFloat(long address, float val) {
  union {
    float val;
    uint8_t byte[4];
  } read;
  const uint8_t* p = (const uint8_t*)(const void*)&val;
  flash.writeBytes(address, p, 4);
  flash.readBytes(address, read.byte, 4);
  return read.val == val;
}

char buffer[20];
int bidx = 0;
long address = -1;
boolean nextInt = false;
boolean nextFloat = false;

void zeroBuf() {
  while (bidx > 0) {
    bidx--;
    buffer[bidx] = ' ';
  }
}

void loop() {
  if (!Serial.available()) { return; }
  char c = Serial.read();
  if (bidx == 0 && c == 'e') {
    flash.blockErase32K(0);
    nextInt = nextFloat = false;
    zeroBuf();
  } else if (bidx > 0 && c == 'i') {
    nextInt = true;
    address = atol(buffer);
    zeroBuf();
  } else if (bidx > 0 && c == 'f') {
    nextFloat = true;
    address = atol(buffer);
    zeroBuf();
  } else if (c == '\n' && nextInt) {
    int val = atoi(buffer);
    nextInt = nextFloat = false;
    zeroBuf();
    if (!writeInt(address, val)) {
      Serial.print("flash write int error, address ");
      Serial.println(address);
      return;
    }
    Serial.print(address);
    Serial.print("i");
    Serial.println(val, DEC);
  } else if (c == '\n' && nextFloat) {
    float val = atof(buffer);
    zeroBuf();
    nextInt = nextFloat = false;
    if (!writeFloat(address, val)) {
      Serial.print("flash write float error, address ");
      Serial.println(address);
      return;
    }
    Serial.print(address);
    Serial.print("f");
    Serial.println(val, 4);
  } else {
    buffer[bidx++] = c;
  }
}

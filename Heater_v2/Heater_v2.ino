#include <Adafruit_MAX31855.h>
#include <SPI.h>

#define heater 13

#define MAXDO1 3
#define MAXCS1 4
#define MAXCLK1 5

#define MAXCS2 7

Adafruit_MAX31855 thermocouple(MAXCLK1, MAXCS1, MAXDO1);
Adafruit_MAX31855 thermocouple2(MAXCLK1, MAXCS2, MAXDO1);

double last_temps[10];
int counter;
double old_dT;

void setup() {
  Serial.begin(9600);
  Serial.println("HEATER INITIALIZE");
  pinMode(heater, OUTPUT);
  digitalWrite(heater, HIGH);
  delay(500);
}

void loop() {
  Serial.print("Ice Temp = ");
  double c = thermocouple.readCelsius();
  if (isnan(c)) {
    Serial.println("Something wrong with thermocouple!");
  } else {
    Serial.println(c);
  }
  Serial.print("Heater Temp = ");
  double c2 = thermocouple2.readCelsius() - 0.5;
  if (isnan(c2)) {
    Serial.println("Something wrong with thermocouple!");
  } else {
    Serial.println(c2);
  }
  double dT = run_avg(c2 - c);
  Serial.print("Delta T = ");
  Serial.println(dT);
  Serial.print("dT/dS = ");
  Serial.println(dT - old_dT);
  Serial.println();
  old_dT = dT;

  delay(1000);
}

double run_avg(double temp) {
  last_temps[counter] = temp;
  double avg = 0;
  for (int i = 0; i < 10; i++) {
    avg += last_temps[i];
  }
  counter = counter + 1;
  if (counter == 10) {
    counter = 0;
  }
  avg = avg / 10;
  return avg;
}

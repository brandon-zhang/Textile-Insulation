#include <Adafruit_MAX31855.h>
#include <SPI.h>

#define heater 13

#define MAXDO1 3
#define MAXCS1 4
#define MAXCLK1 5

#define MAXCS2 7

Adafruit_MAX31855 thermocouple(MAXCLK1, MAXCS1, MAXDO1);
Adafruit_MAX31855 thermocouple2(MAXCLK1, MAXCS2, MAXDO1);

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
   double c2 = thermocouple2.readCelsius();
   if (isnan(c2)) {
     Serial.println("Something wrong with thermocouple!");
   } else {
     Serial.println(c2);
   }
   if (c2 > 100){
    digitalWrite(heater, LOW);
   } else {
    digitalWrite(heater, HIGH);
   }
   Serial.print("Delta T = ");
   Serial.println(c2-c);
   //Serial.print("F = ");
   //Serial.println(thermocouple.readFarenheit());
 
   delay(1000);
}

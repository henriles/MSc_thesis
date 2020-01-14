// Code for controlling sample rotator with arduino using L298 stepper motor driver card
// Henri Leskinen 14.1.2020
// University of Eastern Finland

#include <Keypad.h>
#include <Wire.h>
#include <LiquidCrystal_I2C.h>
#include <Stepper.h>
#define VReadPin 0
#define ENA A5
#define ENB A4

const int StepsPerRev = 200;
const byte rows = 4; //four rows
const byte cols = 4; //four cols
// define keypad layout
char keys[rows][cols] = {
  {'1', '2', '3', 'A'},
  {'4', '5', '6', 'B'},
  {'7', '8', '9', 'C'},
  {'*', '0', '-', 'D'}
};
byte rowPins[rows] = {9, 10, 11, 12}; //row pinouts of the keypad
byte colPins[cols] = {5, 6, 7, 8}; //column pinouts of the keypad

Keypad keypad = Keypad( makeKeymap(keys), rowPins, colPins, rows, cols );
LiquidCrystal_I2C lcd(0x27, 20, 4); // lcd: 4x20 characters and it's hex-adress is 0x27
Stepper stepmot(StepsPerRev, A3, 4, A1, A0); // pin A2 is kaput

int curlox = 0;
int curloy = 0;
byte mode = 0;
byte backto = 0; // 1 = auto, 2 = manual
int keycount = 0;
char charray[4];
int cur_angle = 0;
int step_size = 0;
boolean wait_triggers = false;
int trigger_count = 0;
int ang_array[8];
int ang_ind = 0;
volatile boolean rotation_flag = false;

byte customChar[8] = {
  B01100,
  B10010,
  B10010,
  B01100,
  B00000,
  B00000,
  B00000,
  B00000
};

void setup()
{
  pinMode(VReadPin, INPUT_PULLUP); // define trigger channel
  pinMode(ENA, OUTPUT);
  pinMode(ENB, OUTPUT);
  stepmot.setSpeed(40); // speed = 20 -> 75 sec per 180 deg sample rotation
  // speed 40 -> 37 sec / 180 deg
  Serial.begin(9600); // Start the keypad
  //lcd.init();
  lcd.begin(); // Start thr lcd
  lcd.backlight();
  lcd.createChar(0, customChar);
  // Define interrupt protocol fot trigger channel
  attachInterrupt(digitalPinToInterrupt(VReadPin), InterruptSR, RISING);
  Reset();

}

void loop()
{
  // if trigger was received, turn the sample.
  if (rotation_flag) {
    Next_Angle();
    rotation_flag = false;
  }

  char key = keypad.getKey();

  if (backto == 3) {
    key = 'A';
  } else if (backto == 2 || backto == 4) {
    key = 'B';
  }

  // if a key was pressed in the keypad, go trough this logic to decide what to do.
  if (key != NO_KEY) {

    switch (key) {
      case 'A':
        if (mode == 0 ) {
          // if in main menu, open auto rotation -menu
          mode = 1; // auto-mode
          lcd.clear();
          lcdPrintLoc(1, 0, "A: Const step size");
          lcdPrintLoc(2, 0, "B: Define angles");
        } else if (mode == 1 || backto == 3) {
          // if in auto rotation -menu, choose const step size
          mode = 3;
          backto = 0; // not returning by default
          lcd.clear();
          lcd.print("AUTO: Set step size");
          lcdPrintLoc(1, 0, "(deg) and press C");
          lcd.setCursor(5, 3);
          lcd.write(byte(0));
          lcd.setCursor(0, 3);
          lcd.blink_on();
        }
        break;





      case 'B':
        if (mode == 0 || backto == 2) {
          // if coming from main menu or coming from manual rotation, go to manual rotation
          mode = 2; // manual mode
          backto = 0; // not returning to main menu by default
          lcd.clear();
          lcd.print("MANUAL     D = EXIT");
          lcdPrintLoc(1, 0, "Set angle & press C");
          lcdPrintLoc(2, 12, "Cur:");
          lcdPrintLoc(3, 12, String(String(cur_angle) ));
          lcd.write(byte(0));
          lcd.setCursor(5, 3);
          lcd.write(byte(0));
          lcd.setCursor(0, 3);
          lcd.blink_on();
        } else if (mode == 1 || backto == 4) {
          // if coming from auto rotation -menu, choose Auto rotation and "define steps"
          mode = 4;
          backto = 0;
          lcd.clear();
          lcdPrintLoc(0, 0, "Insert steps");
          lcdPrintLoc(1, 0, "* = OK, C = Run");
          lcd.setCursor(5, 3);
          lcd.write(byte(0));
          lcd.setCursor(0, 3);
          lcd.blink_on();
        }
        break;




      case 'C':
        // 'C' means OK
        // cant press "OK" if nothing is inserted
        if ((mode == 1 || mode == 3) && keycount == 0) {
        }
        else if (mode == 0) {
          // if in main menu, define current angle to new 0.
          cur_angle = 0;
          lcd.clear();
          Reset();
        } else if (mode == 3) {
          // if in "constant step size" menu. Save the inserted step size
          // and start waiting triggers
          lcd.blink_off();
          lcd.clear();
          lcd.print(String("AUTO " + String(trigger_count)) + " trgs, " + String(cur_angle));
          lcd.write(byte(0));
          step_size = CalAngle(0);
          if (step_size == 0) {
            lcdPrintLoc(1, 0, String("Invalid angle "));
            backto = 3;
          } else {
            lcdPrintLoc(1, 0, String("Step size " + String(step_size) ));
            lcd.write(byte(0));
            lcdPrintLoc(2, 0, "Waiting for triggers");
            lcdPrintLoc(3, 0, "Press D to quit");
            wait_triggers = true;
          }
          delay(2000);
        } else if (mode == 4 && wait_triggers == false) {
          // if in "define angles" -menu, start waiting for triggers
          lcd.blink_off();
          lcd.clear();
          lcdPrintLoc(0, 0, String("AUTO 0 trgs, " + String(cur_angle)));
          lcd.write(byte(0));
          lcdPrintLoc(1, 0, "Waiting trigs,D=Quit");
          lcd.setCursor(0, 2);
          for (int i = 0; i < ang_ind; i++) {
            if (i == 4) {
              lcd.setCursor(0, 3);
            }
            lcd.print(String(ang_array[i]) + " ");
          }
          wait_triggers = true;
        } else if (mode == 2) {
          // if in manual, rotate the stepper (sample)
          lcd.blink_off();
          lcd.clear();
          lcd.print("MANUAL     D = EXIT");
          int new_angle = CalAngle(cur_angle);
          if (new_angle == cur_angle) {
            lcdPrintLoc(2, 0, String("Invalid angle "));
            delay(3000);
          } else {
            lcdPrintLoc(2, 0, String("Rotating to " + String(new_angle)));
            RotateTo(new_angle);
            lcd.write(byte(0));
            cur_angle = new_angle;
          }
          backto = 2;
        }
        keycount = 0;
        break;

      case 'D':
        // Go to main menu if not already there
        Serial.println(String("keycount - mode - wait_triggers:" + String(keycount)
        + String(mode) + String(wait_triggers)));
        if (mode != 0) {
          lcd.clear();
          Reset();
        }
        break;

      case '-':
        if (mode > 1 && keycount == 0) {
          lcd.print(key);
          charray[keycount] = key;
          keycount++;
        } else if (mode > 1) {
          Erase();
        }
        break;

      case '*':
        // if inserting angles manually in 4th mode, add it.
        Serial.println(String("keycount - mode - wait_triggers:" + String(keycount)
        + String(mode) + String(wait_triggers))); // this is just for debuggind purposes
        if (mode == 4 && keycount != 0 && ang_ind < 8) {
          ang_array[ang_ind] = CalAngle(0);
          if (ang_ind < 5) {
            lcdPrintLoc(2, ang_ind * 4, String(ang_array[ang_ind]));
          } else if (ang_ind > 4) {
            lcdPrintLoc(3, 6 + ((ang_ind - 5) * 4), String(ang_array[ang_ind]));
          }
          Erase();
          ang_ind++;
          Serial.println(String( "ang_array(i) / cal angle:" + String(ang_array[ang_ind - 1]) ));
        }
        break;

      default:
        // if a number was pressed in correct menu, print it.
        if (mode > 1 && keycount < 4 && !wait_triggers) {
          lcd.print(key);
          charray[keycount] = key;
          keycount++;
          Serial.println(String(String(key) + " keycount: " + String(keycount) ));
        }
        break;
    }
  }
}

int ch2int(char c) {
  int db = (c - '0');
  return db;
}

int CalAngle(int old_ang) {
  // transform inserted charactes from character array to integer
  int val = 0;
  boolean neg = false;
  // if first character is '-', the value is negative
  if (charray[0] == '-') {
    Serial.println(String("Char-array: " + String(charray) ));
    charray[0] = charray[1];
    charray[1] = charray[2];
    charray[2] = charray[3];
    keycount--;
    neg = true;
  }

  int temp = 0;
  // pow(a,b) doesnt work properly for integers so i had to write my own
  for (int i = 0; i < keycount; i++) {
    int powi = 1;
    for (int j = 0; j < keycount - (i + 1); j++) {
      powi = powi * 10;
    }
    val = val + ch2int(charray[i]) * powi;
  }


  if (neg) {
    val = (-1) * val;
  }

  if (abs(val) > 190) {
    val = old_ang;
  }

  //int value = int(val);
  Serial.println(String("angle:" + String(val)));
  return val;
}

void Reset () {
  mode = 0;
  backto = 0;
  keycount = 0;
  trigger_count = 0;
  ang_ind = 0;
  wait_triggers = false;
  lcd.blink_off();
  lcdPrintLoc(0, 0, String("MAIN      Cur: " + String(cur_angle) ));
  lcd.write(byte(0));
  lcdPrintLoc(1, 0, "A: Auto rotation");
  lcdPrintLoc(2, 0, "B: Manual rotation");
  lcdPrintLoc(3, 0, "C: Set this angle 0");
  lcd.write(byte(0));
  digitalWrite(ENA, LOW);
  digitalWrite(ENB, LOW);
}

void lcdPrintLoc (int row, int col, String st) {
  lcd.setCursor(col, row);
  lcd.print(st);
}

void RotateTo (int th) {
  // wake up the stepper and rotete it
  digitalWrite(ENA, HIGH);
  digitalWrite(ENB, HIGH);
  double stepsPerDeg = 29;
  int steps = (double(th - cur_angle)) *  stepsPerDeg;
  //steps = double(th - cur_angle);
  stepmot.step(steps);
  delay(1000);
  digitalWrite(ENA, LOW);
  digitalWrite(ENB, LOW);
}

// if trigger is received, turn on a flag
void InterruptSR () {
  if (wait_triggers) {
    rotation_flag = true;
  }
}

void Next_Angle() {
  if (mode == 3) {
    //Serial.println(digitalRead(VReadPin));
    if (abs(cur_angle + step_size) <= 180) {
      delay(abs(10 * step_size));
      RotateTo(cur_angle + step_size);
      cur_angle = cur_angle + step_size;
      trigger_count++;
      lcdPrintLoc(0, 0, "                  ");
      lcdPrintLoc(0, 0, String("AUTO " + String(trigger_count)) + " trgs, "
      + String(cur_angle));
      lcd.write(byte(0));
      delay(10);

    }
  } else if (mode == 4) {
    //Serial.println(digitalRead(VReadPin));
    if (trigger_count < ang_ind) {
      delay(abs(10 * ang_array[trigger_count]));
      RotateTo(ang_array[trigger_count]);
      cur_angle = ang_array[trigger_count];
      trigger_count++;
      lcdPrintLoc(0, 0, "                  ");
      lcdPrintLoc(0, 0, String("AUTO " + String(trigger_count)) + " trgs, "
      + String(cur_angle));
      lcd.write(byte(0));
      delay(10);
    }

  }
}

void Erase () {
  char charray[4];
  lcdPrintLoc(3, 0, "     ");
  lcd.setCursor(0, 3);
  keycount = 0;
}

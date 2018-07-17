BufferedReader reader1, reader2, reader3;
String lineX, lineY, lineZ;
String sciezka;

int xsize=1000, ysize=1000, counter, choice, anim, liczba_czastek, liczba_warstw, szerokosc, wskA, wskC;
float skala=200;
float x0, y0, z0, eyeX, eyeY, eyeZ;
float teta, fi, r, deltaAlfa, deltaR, sphR;
float teta0, fi0, r0, sphR0;
float rotX, rotY, rotZ;

void setup() {
  counter=0;
  anim=1;
  choice=2;
  
  liczba_warstw = 12;
  reader1 = createReader("RX.txt"); 
  reader2 = createReader("RY.txt"); 
  reader3 = createReader("RZ.txt"); 
  frameRate(15);
  size(600, 600, P3D);
  background(255);  
  r0=200;
  teta0=0;
  fi0=0;
  sphR0=5;

  if (choice==1)
  {
    teta0=6.28*0.25;
    r0=200;
  }

  r=r0;
  teta=teta0;
  fi=fi0;
  sphR=sphR0;

  x0=0.5*xsize;
  y0=0.5*ysize;
  z0=0.5*xsize;
  deltaAlfa=0.05;
  deltaR=50;
}
void draw() 
{
  counter=counter+1;

  try 
  {
    lineX = reader1.readLine();
    lineY = reader2.readLine();
    lineZ = reader3.readLine();
  } 
  catch (IOException e) 
  {
    e.printStackTrace();
    lineX = null;
    lineY = null;
    lineZ = null;
  }
  if (lineX == null || lineY == null || lineZ == null) 
  {
    noLoop();
  } else 
  {
    String[] piecesX = split(lineX, " ");
    String[] piecesY = split(lineY, " ");
    String[] piecesZ = split(lineZ, " ");
    liczba_czastek=piecesX.length;
    szerokosc=int(sqrt(float(liczba_czastek)/(liczba_warstw*4)));
    wskA=szerokosc*szerokosc*4*4;
    wskC=liczba_czastek-wskA;

    float[] x= new float[liczba_czastek];
    float[] y= new float[liczba_czastek];
    float[] z= new float[liczba_czastek];


    for (int i=0; i<liczba_czastek; i++)
    {
      x[i] = float(piecesX[i])*20+500;
      y[i] = float(piecesY[i])*20+500;
      z[i] = float(piecesZ[i])*20+500;
    }


    background(255);  
    noStroke();
    lights();

    if (keyPressed) 
    {
      if (key == 'a') 
      {
        teta=teta-deltaAlfa;
      } else if (key == 'd') 
      {
        teta=teta+deltaAlfa;
      } else if (key == 'w') 
      {
        fi=fi-deltaAlfa;
      } else if (key == 's') 
      {
        fi=fi+deltaAlfa;
      } else if (key == 'q') 
      {
        r=r+deltaR;
      } else if (key == 'e') 
      {
        r=r-deltaR;
      } else if (key == 'r') 
      {
        sphR=sphR+0.5;
      } else if (key == 'f') 
      {
        sphR=sphR-0.5;
      } else if (key == 'c') 
      {
        r=r0;
        teta=teta0;
        fi=fi0;
        sphR=sphR0;
        rotX=0;
        rotY=0;
        rotZ=0;
      } else if (key == 't') 
      {
        rotY=rotY+0.25;
      } else if (key == 'y') 
      {
        rotY=rotY-0.25;
      }
    } 

    eyeX=x0+r*cos(teta)*cos(fi);
    eyeY=y0+r*cos(teta)*sin(fi);
    eyeZ=z0+r*sin(teta);


    if (choice==2)
    {
      camera(eyeX, eyeY, eyeZ, // eyeX, eyeY, eyeZ
        0.5*xsize, 0.5*ysize, 0.5*ysize, // centerX, centerY, centerZ
        0.0, 1.0, 0.0); // upX, upY, upZ
    }

    if (choice==3)
    {
      ortho(-0.2*xsize, 0.6*xsize, -0.6*ysize, 0.2*ysize);    
      rotateX(rotX);
      rotateY(rotY);
      rotateZ(rotZ);
    }

    if (choice==1)
    {

      fill(178, 34, 34);
      for (int i=0; i<wskA; i++)
      {
        pushMatrix();
        translate(x[i], z[i], y[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }

      fill(50, 205, 50);
      for (int i=wskA; i<(wskC); i++)
      {
        pushMatrix();
        translate(x[i], z[i], y[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }

      fill(178, 34, 34);
      for (int i=(wskC); i<liczba_czastek; i++)
      {
        pushMatrix();
        translate(x[i], z[i], y[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }
    } else if (choice==2)
    {
      teta=6.28/4;
      fill(255, 0, 0);
      for (int i=0; i<wskA; i++)
      {
        ellipse(y[i], z[i], sphR, sphR);
      }

      fill(0, 255, 0);
      for (int i=wskA; i<wskC; i++)
      {
        ellipse(y[i], z[i], sphR, sphR);
      }

      fill(255, 0, 0);
      for (int i=(wskC); i<liczba_czastek; i++)
      {
        ellipse(y[i], z[i], sphR, sphR);
      }

      if (counter>2000 && anim==1)
      {
        println(counter);
      }
    } else if (choice==3)
    {

      fill(178, 34, 34);
      for (int i=0; i<wskA; i++)
      {
        pushMatrix();
        translate(y[i], z[i], x[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }

      fill(50, 205, 50);
      for (int i=wskA; i<wskC; i++)
      {
        pushMatrix();
        translate(y[i], z[i], x[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }

      fill(178, 34, 34);
      for (int i=(wskC); i<liczba_czastek; i++)
      {
        pushMatrix();
        translate(y[i], z[i], x[i]);
        sphereDetail(10);
        sphere(sphR);
        popMatrix();
      }
      println(counter);
    }
  }
} 
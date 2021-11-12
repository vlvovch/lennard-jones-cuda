#include "MDSystemGL.h"

//#include <GL/glew.h>
//#include <gl\gl.h>
//#include <qopenglextensions_p.h>
#include <QDebug>
#include <QGLWidget>
//#include <GL/glew.h>
//#include <gl\glut.h>
//#include <gl\glu.h>

MDSystemGL::MDSystemGL(const MDSystem::MDSystemConfiguration& config, bool renderSpheres) : m_renderSpheres(renderSpheres)
{
  isDList = 0;
  m_system = new MDSystem(config);

  usepbo = false;

  camera_trans[0] = camera_trans[1] = camera_trans[2] = 0;
  camera_rot[0] = camera_rot[1] = camera_rot[2] = 0;
  camera_trans_lag[0] = camera_trans_lag[1] = camera_trans_lag[2] = 0;
  camera_rot_lag[0] = camera_rot_lag[1] = camera_rot_lag[2] = 0;
}

void MDSystemGL::Reinitialize(const MDSystem::MDSystemConfiguration& config, bool renderSpheres)
{
  m_renderSpheres = renderSpheres;
  isDList = 0;
  m_system->Reinitialize(config);
}

void MDSystemGL::draw()
{
  //qDebug() << iter++ << "\n";
  if (!isDList) initSphere();
  glColor3ub(255, 255, 255);
  glColor3ub(0, 0, 0);
  glPushMatrix();
  glTranslatef(0, 0, -300.0);

  for (int i = 0; i < 3; ++i)
  {
    camera_trans_lag[i] += (camera_trans[i] - camera_trans_lag[i]);// * inertia;
    camera_rot_lag[i] += (camera_rot[i] - camera_rot_lag[i]);// * inertia;
  }
  glTranslatef(camera_trans_lag[0],
    camera_trans_lag[1],
    camera_trans_lag[2]);
  //glRotatef(-90.0, 1.0, 0.0, 0.0);
  glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
  glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);

  /*glutWireCube(150.0);
  glTranslatef(-75.0, -75.0, -75.0);*/
  //glTranslatef(-75.0, -75.0, -75.0);
  glEnable(GL_LIGHTING);
  glTranslatef(-75.0, -75.0, -75.0);
  glLineWidth(2.0f);
  glBegin(GL_LINES);
  glVertex3f(0., 0., 0.);
  glVertex3f(150., 0., 0.);

  glVertex3f(0., 0., 0.);
  glVertex3f(0., 150., 0.);

  glVertex3f(150., 0., 0.);
  glVertex3f(150., 150., 0.);

  glVertex3f(0., 150., 0.);
  glVertex3f(150., 150., 0.);

  glVertex3f(0., 0., 150.);
  glVertex3f(150., 0., 150.);

  glVertex3f(0., 0., 150.);
  glVertex3f(0., 150., 150.);

  glVertex3f(150., 0., 150.);
  glVertex3f(150., 150., 150.);

  glVertex3f(0., 150., 150.);
  glVertex3f(150., 150., 150.);

  glVertex3f(0., 0., 0.);
  glVertex3f(0., 0., 150.);

  glVertex3f(150., 0., 0.);
  glVertex3f(150., 0., 150.);

  glVertex3f(0., 150., 0.);
  glVertex3f(0., 150., 150.);

  glVertex3f(150., 150., 0.);
  glVertex3f(150., 150., 150.);


  glEnd();

  double fraction = 0.4;

  //subcube
  glPushAttrib(GL_ENABLE_BIT);
  glLineWidth(3.0f);
  glLineStipple(1, 8000);
  glColor3ub(255, 0, 0);
  glEnable(GL_LINE_STIPPLE);
  glBegin(GL_LINES);
  for (double dir = -0.5; dir < 0.6; dir += 1.) {
    glVertex3f(0., 0., 75. + dir * 150. * fraction);
    glVertex3f(150., 0., 75. + dir * 150. * fraction);

    glVertex3f(0., 150., 75. + dir * 150. * fraction);
    glVertex3f(150., 150., 75. + dir * 150. * fraction);

    glVertex3f(0., 0., 75. + dir * 150. * fraction);
    glVertex3f(0., 150., 75. + dir * 150. * fraction);

    glVertex3f(150., 0., 75. + dir * 150. * fraction);
    glVertex3f(150., 150., 75. + dir * 150. * fraction);
  }
  glEnd();
  glPopAttrib();

  //glTranslatef(75.0, 75.0, 75.0);
  glColor3ub(255, 127, 0);
  //glColor3ub(0, 0, 0);
  /*glBegin( GL_POINTS );
  for(int i=0;i<N;++i)
  {
      glVertex3f( x[i] * 150.0 / sz, y[i] * 150.0 / sz, z[i] * 150.0 / sz );
  }
  glEnd();*/
  if (m_renderSpheres)
  {
    for (int i = 0; i < 4 * m_system->m_config.N; i += 4)
    {
      double z = m_system->h_Pos[i + 2];
      if (std::abs(z - 0.5 * m_system->L) < 0.5 * fraction * m_system->L) {
        //glColor3ub(0, 0, 255);
        glColor3ub(255, 0, 0);
      }
      else {
        //glColor3ub(255, 0, 0);
        glColor3ub(30, 30, 30);
      }
      
      glPushMatrix();
      glTranslatef(m_system->h_Pos[i] * 150.0 / m_system->L, m_system->h_Pos[i + 1] * 150.0 / m_system->L, m_system->h_Pos[i + 2] * 150.0 / m_system->L);
      //glScalef(0.5, 0.5, 0.5);
      glCallList(dListSphere);
      //glutSolidSphere( 0.5 * 150.0 / sz, 20, 20);
      glPopMatrix();
    }
    glDisable(GL_LIGHTING);
  }
  else
  {
    glDisable(GL_LIGHTING);
    glEnable(GL_POINT_SPRITE_ARB);
    glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
    glPointSize(2.0f);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
    //glDepthMask(GL_TRUE);
    //glEnable(GL_DEPTH_TEST);

    /*glUseProgram(m_program);
    GLuint texLoc = glGetUniformLocation(m_program, "splatTexture");
    glUniform1i(texLoc, 0);

    glActiveTextureARB(GL_TEXTURE0_ARB);
    glBindTexture(GL_TEXTURE_2D, m_texture);*/

    glUseProgram(m_program);
    GLuint texLoc = glGetUniformLocation(m_program, "splatTexture");
    glUniform1i(texLoc, 0);
    //glUniform1f( glGetUniformLocation(m_program2, "pointScale"), h / tanf(60.0*0.5f*(float)PI/180.0f) );
    //glUniform1f( glGetUniformLocation(m_program2, "pointRadius"), .5 );
    //glColor3f(1, 1, 1);

    glActiveTextureARB(GL_TEXTURE0_ARB);
    //glActiveTexture(GL_TEXTURE0_ARB);
    glBindTexture(GL_TEXTURE_2D, m_texture);

    glEnableClientState(GL_VERTEX_ARRAY);


    glVertexPointer(4, GL_FLOAT, 0, m_system->h_Pos);
    glDrawArrays(GL_POINTS, 0, m_system->m_config.N);
    glDisableClientState(GL_VERTEX_ARRAY);

    glUseProgram(0);

    //glDisable(GL_POINT_SPRITE_ARB);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
  }
  glPopMatrix();
}


const char vertexShader[] =
{
    "void main()                                                            \n"
    "{                                                                      \n"
    "    float pointSize = 500.0 * gl_Point.size;                           \n"
    "    vec4 vert = gl_Vertex;												\n"
  //"    vert.w = 1.0;														\n"
  "    vec3 pos_eye = vec3 (gl_ModelViewMatrix * vert);                   \n"
  "    gl_PointSize = max(1.0, pointSize / (1.0 - pos_eye.z));            \n"
  "    gl_TexCoord[0] = gl_MultiTexCoord0;                                \n"
  //"    gl_TexCoord[1] = gl_MultiTexCoord1;                                \n"
  "    gl_Position = ftransform();                                        \n"
  //"gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;  \n"
  "    gl_FrontColor = gl_Color;                                          \n"
  "}                                                                      \n"
};

const char pixelShader[] =
{
    "uniform sampler2D splatTexture;                                        \n"

    "void main()                                                            \n"
    "{                                                                      \n"
    "    vec4 color = (0.6 + 0.4 * gl_Color) * texture2D(splatTexture, gl_TexCoord[0].st); \n"
    "    gl_FragColor =                                                     \n"
    "         color * lerp(vec4(0.1, 0.0, 0.0, color.w), vec4(1.0, 0.7, 0.3, color.w), color.w);\n"
    "}                                                                      \n"
};

extern const char* vertexShader2;
extern const char* spherePixelShader2;

void MDSystemGL::InitGL()
{
  //glewInit();
  QOpenGLFunctions::initializeOpenGLFunctions();
  QOpenGLExtension_ARB_multitexture::initializeOpenGLFunctions();

  m_vertexShader = glCreateShader(GL_VERTEX_SHADER);
  m_pixelShader = glCreateShader(GL_FRAGMENT_SHADER);

  const char* v = vertexShader;
  const char* p = pixelShader;
  glShaderSource(m_vertexShader, 1, &v, 0);
  glShaderSource(m_pixelShader, 1, &p, 0);

  glCompileShader(m_vertexShader);
  glCompileShader(m_pixelShader);

  m_program = glCreateProgram();

  glAttachShader(m_program, m_vertexShader);
  glAttachShader(m_program, m_pixelShader);

  glLinkProgram(m_program);

  _createTexture(32);

  GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
  GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

  glShaderSource(vertexShader, 1, &vertexShader2, 0);
  glShaderSource(fragmentShader, 1, &spherePixelShader2, 0);

  glCompileShader(vertexShader);
  glCompileShader(fragmentShader);

  m_program2 = glCreateProgram();

  glAttachShader(m_program2, vertexShader);
  glAttachShader(m_program2, fragmentShader);

  glLinkProgram(m_program2);

  glEnable(GL_DEPTH_TEST);

  glEnable(GL_COLOR_MATERIAL);

  GLfloat afPropertiesAmbient[] = { 0.50, 0.50, 0.50, 1.00 };
  GLfloat afPropertiesDiffuse[] = { 0.75, 0.75, 0.75, 1.00 };
  GLfloat afPropertiesSpecular[] = { 1.00, 1.00, 1.00, 1.00 };
  //GLfloat afSpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };
  GLfloat afSpecularWhite[] = { 0.50, 0.50, 0.50, 1.00 };

  glEnable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glLightfv(GL_LIGHT0, GL_AMBIENT, afPropertiesAmbient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, afPropertiesDiffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, afPropertiesSpecular);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

  /*float pos[4] = {20,-20,-200,1};
  float dir[3] = {1,1,-1};
  float color[4] = {1,1,1,1};
  float mat_specular[] = {1,1,1,1};

  glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
  glLightfv(GL_LIGHT0, GL_POSITION, pos);
  glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, dir);*/
  //float dir[3] = {1,0,0};
  /*float pos[4] = {0,-1,0,0};
  glLightfv(GL_LIGHT0, GL_POSITION, pos);*/

  glEnable(GL_LIGHT0);

  //glMaterialfv(GL_BACK,  GL_AMBIENT,   afAmbientGreen);
  //glMaterialfv(GL_BACK,  GL_DIFFUSE,   afDiffuseGreen);
  //glMaterialfv(GL_FRONT, GL_AMBIENT,   afAmbientRed);
  //glMaterialfv(GL_FRONT, GL_DIFFUSE,   afDiffuseRed);
  //glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
  glMaterialfv(GL_FRONT, GL_SPECULAR, afSpecularWhite);
  glMaterialf(GL_FRONT, GL_SHININESS, 40.0);
  //glMaterialf( GL_FRONT, GL_SHININESS, 4000.0);

  glClearColor(0, 0, 0, 1);
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPointSize(2.0);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glutPostRedisplay();
}


void MDSystemGL::initSphere()
{
  dListSphere = glGenLists(1);
  glNewList(dListSphere, GL_COMPILE);
  //GLUquadricObj *quadObj;
  //quadObj = gluNewQuadric();
  //gluSphere(quadObj, 0.5 * 150.0 / m_system->L, 20, 20);
  //gluDeleteQuadric( quadObj );
  float PI = 3.14159265358;
  int iterPhi = 20;
  int iterTheta = 20;
  float dphi = (float)(2.f * PI / iterPhi);
  float dtheta = (float)(PI / iterTheta);
  float phi = 0.f, theta = 0.f;
  float rad = 0.5 * 150.0 / m_system->L * 0.25;
  float x[4];
  float y[4];
  float z[4];
  glBegin(GL_TRIANGLES);
  for (int i = 0; i < iterPhi; ++i) {
    float phi = (float)(PI / 2.f + i * dphi);
    for (int j = 0; j < iterTheta; ++j) {
      theta = j * dtheta;
      x[0] = (float)cos(phi) * (float)sin(theta);
      y[0] = (float)sin(phi) * (float)sin(theta);
      z[0] = (float)cos(theta);
      x[1] = (float)cos(phi) * (float)sin(theta + dtheta);
      y[1] = (float)sin(phi) * (float)sin(theta + dtheta);
      z[1] = (float)cos(theta + dtheta);
      x[2] = (float)cos(phi + dphi)
        * (float)sin(theta + dtheta);
      y[2] = (float)sin(phi + dphi)
        * (float)sin(theta + dtheta);
      z[2] = (float)cos(theta + dtheta);
      x[3] = (float)cos(phi + dphi) * (float)sin(theta);
      y[3] = (float)sin(phi + dphi) * (float)sin(theta);
      z[3] = (float)cos(theta);

      glNormal3f(x[0], y[0], z[0]);
      glVertex3f(rad * x[0], rad * y[0], rad * z[0]);
      glNormal3f(x[1], y[1], z[1]);
      glVertex3f(rad * x[1], rad * y[1], rad * z[1]);
      glNormal3f(x[2], y[2], z[2]);
      glVertex3f(rad * x[2], rad * y[2], rad * z[2]);

      glNormal3f(x[0], y[0], z[0]);
      glVertex3f(rad * x[0], rad * y[0], rad * z[0]);
      glNormal3f(x[2], y[2], z[2]);
      glVertex3f(rad * x[2], rad * y[2], rad * z[2]);
      glNormal3f(x[3], y[3], z[3]);
      glVertex3f(rad * x[3], rad * y[3], rad * z[3]);
    }
  }
  glEnd();

  glEndList();
}


//------------------------------------------------------------------------------
// Function     	  : EvalHermite
// Description	    :
//------------------------------------------------------------------------------
/**
* EvalHermite(float pA, float pB, float vA, float vB, float u)
* @brief Evaluates Hermite basis functions for the specified coefficients.
*/
inline float evalHermite(float pA, float pB, float vA, float vB, float u)
{
  float u2 = (u * u), u3 = u2 * u;
  float B0 = 2 * u3 - 3 * u2 + 1;
  float B1 = -2 * u3 + 3 * u2;
  float B2 = u3 - 2 * u2 + u;
  float B3 = u3 - u;
  return(B0 * pA + B1 * pB + B2 * vA + B3 * vB);
}


unsigned char* createGaussianMap(int N)
{
  float* M = new float[2 * N * N];
  unsigned char* B = new unsigned char[4 * N * N];
  float X, Y, Y2, Dist;
  float Incr = 2.0f / N;
  int i = 0;
  int j = 0;
  Y = -1.0f;
  //float mmax = 0;
  for (int y = 0; y < N; y++, Y += Incr)
  {
    Y2 = Y * Y;
    X = -1.0f;
    for (int x = 0; x < N; x++, X += Incr, i += 2, j += 4)
    {
      Dist = (float)sqrtf(X * X + Y2);
      if (Dist > 1) Dist = 1;
      M[i + 1] = M[i] = evalHermite(1.0f, 0, 0, 0, Dist);
      B[j + 3] = B[j + 2] = B[j + 1] = B[j] = (unsigned char)(M[i] * 255);
    }
  }
  delete[] M;
  return(B);
}

void MDSystemGL::_createTexture(int resolution)
{
  unsigned char* data = createGaussianMap(resolution);
  glGenTextures(1, &m_texture);
  glBindTexture(GL_TEXTURE_2D, m_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, resolution, resolution, 0,
    GL_RGBA, GL_UNSIGNED_BYTE, data);

}

#define STRINGIFY(A) #A

// vertex shader
const char* vertexShader2 = STRINGIFY(
  uniform float pointRadius;  // point size in world space
uniform float pointScale;   // scale to calculate size in pixels
uniform float densityScale;
uniform float densityOffset;
void main()
{
  // calculate window-space point size
  vec3 posEye = vec3(gl_ModelViewMatrix * gl_Vertex);
  float dist = length(posEye);
  gl_PointSize = pointRadius * (pointScale / dist);

  gl_TexCoord[0] = gl_MultiTexCoord0;
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

  gl_FrontColor = gl_Color;
}
);

// pixel shader for rendering points as shaded spheres
const char* spherePixelShader2 = STRINGIFY(
  void main()
{
  const vec3 lightDir = vec3(0.577, 0.577, 0.577);

  // calculate normal from texture coordinates
  vec3 N;
  N.xy = gl_TexCoord[0].xy * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
  float mag = dot(N.xy, N.xy);
  if (mag > 1.0) discard;   // kill pixels outside circle
  N.z = sqrt(1.0 - mag);

  // calculate lighting
  float diffuse = max(0.0, dot(lightDir, N));

  gl_FragColor = gl_Color * diffuse;
}
);

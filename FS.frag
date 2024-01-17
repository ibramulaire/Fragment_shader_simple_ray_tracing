#version 450

#define MAX_DIST 100


out vec4 FragColor;
uniform float Largeur;
uniform float Hauteur;

uniform vec3 cameraPosition;

uniform mat4 MODEL;
uniform mat4 MVP;
uniform mat4 VIEW;
uniform mat4 PERSPECTIVE;

mat4 INVMODEL=inverse(MODEL);

vec2 resolution = vec2(Largeur, Hauteur);
vec2 UVtemps= (gl_FragCoord.xy/resolution)-0.5;
float Aspect_ratio=Hauteur/Largeur;
vec2 uv=vec2(UVtemps.x/Aspect_ratio,UVtemps.y);






vec3 hitpoint=vec3(0,0,0);




vec3 normal;

mat4 rotation(vec3 axis, float angle) {
  angle=radians(angle);
    float cosTheta = cos(angle);
    float sinTheta = sin(angle);
    float oneMinusCosTheta = 1.0 - cosTheta;

   
    axis = normalize(axis);

    mat3 rotationMatrix3 = mat3(
        cosTheta + axis.x * axis.x * oneMinusCosTheta,
        axis.x * axis.y * oneMinusCosTheta - axis.z * sinTheta,
        axis.x * axis.z * oneMinusCosTheta + axis.y * sinTheta,

        axis.y * axis.x * oneMinusCosTheta + axis.z * sinTheta,
        cosTheta + axis.y * axis.y * oneMinusCosTheta,
        axis.y * axis.z * oneMinusCosTheta - axis.x * sinTheta,

        axis.z * axis.x * oneMinusCosTheta - axis.y * sinTheta,
        axis.z * axis.y * oneMinusCosTheta + axis.x * sinTheta,
        cosTheta + axis.z * axis.z * oneMinusCosTheta
    );


    mat4 rotationMatrix = mat4(
        rotationMatrix3[0], 0.0,
        rotationMatrix3[1], 0.0,
        rotationMatrix3[2], 0.0,
        vec4(0.0, 0.0, 0.0, 1.0)
    );

    return inverse(rotationMatrix);
}


mat4 translate(vec3 translation) {
  return inverse( mat4(1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              translation.x, translation.y, translation.z, 1.0));
}



float dot2( in vec3 v ) { return dot(v,v); }

float iplane(vec3 pos,vec3 axe,float degree,vec3 position,in vec3 Pnormal,out vec3 ro,out vec3 rd) 
{
   mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
     ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;


     rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);

    float denominator = dot(Pnormal, rd);

    if (abs(denominator) > 0.0001) {
       return (dot(Pnormal, position-ro)) / denominator;
        }
    
    return MAX_DIST;

}


vec2 isphere(vec3 pos,vec3 axe,float degree,in float radius,out vec3 ro,out vec3 rd )
{
    mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;
    rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    float a =1;
  
    float b = 2.0 * dot(ro, rd);

    float c = dot(ro,ro) - radius;
    float delta=(b*b)-4*a*c;

    if(delta<0)
    return vec2(MAX_DIST);

    return vec2((-b - sqrt(delta)) / (2.0 * a),(-b + sqrt(delta)) / (2.0 * a));

}


vec4 icylinder(vec3 pos,vec3 axe,float degree, in vec3 a, in vec3 b, float radiuis,out vec3 ro,out vec3 rd)
{

    mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;
    rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    vec3  ba = b  - a;
    vec3  oc = ro - a;
    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoc = dot(ba,oc);
    float k2 = baba            - bard*bard;
    float k1 = baba*dot(oc,rd) - baoc*bard;
    float k0 = baba*dot(oc,oc) - baoc*baoc - radiuis*radiuis*baba;
    float h = k1*k1 - k2*k0;
    if( h<0.0 ) return vec4(MAX_DIST);//no intersection
    h = sqrt(h);
    float t = (-k1-h)/k2;
    // body
    float y = baoc + t*bard;
    if( y>0.0 && y<baba ) return vec4( t, normalize((oc+t*rd - ba*y/baba)/radiuis ));
    // caps
    t = ( ((y<0.0) ? 0.0 : baba) - baoc)/bard;
    if( abs(k1+k2*t)<h )
    {
        return vec4( t,normalize(ba*sign(y)/sqrt(baba) ));
    }
    return vec4(MAX_DIST);//no intersection
}


vec2 iInfinit_Cylinder( vec3 pos,vec3 axer,float degree, in vec3 base, in vec3 axe, float radius ,out vec3 ro,out vec3 rd)
{   mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axer,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;
    rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    vec3  oc = ro - base;
    float card = dot(axe,rd);
    float caoc = dot(axe,oc);
    float a = 1.0 - card*card;
    float b = dot( oc, rd) - caoc*card;
    float c = dot( oc, oc) - caoc*caoc - radius*radius;
    float h = b*b - a*c;
    if( h<0.0 ) return vec2(MAX_DIST);
    h = sqrt(h);
    return vec2(-b-h,-b+h)/a;
}


// cone defined by extremes pa and pb, and radious ra and rb

vec4 coneIntersect( vec3 pos,vec3 axe,float degree, in vec3 pa, in vec3 pb, in float ra, in float rb,out vec3 ro,out vec3 rd )
{
    mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;
    rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    vec3  ob = ro - pb;
    float m0 = dot(ba,ba);
    float m1 = dot(oa,ba);
    float m2 = dot(rd,ba);
    float m3 = dot(rd,oa);
    float m5 = dot(oa,oa);
    float m9 = dot(ob,ba); 
    
    // caps
    if( m1<0.0 )
    {
        if( dot2(oa*m2-rd*m1)<(ra*ra*m2*m2) ) // delayed division
            return vec4(-m1/m2,normalize(-ba*inversesqrt(m0)));
    }
    else if( m9>0.0 )
    {
    	float t = -m9/m2;                     // NOT delayed division
        if( dot2(ob+rd*t)<(rb*rb) )
            return vec4(t,normalize(ba*inversesqrt(m0)));
    }
    
    // body
    float rr = ra - rb;
    float hy = m0 + rr*rr;
    float k2 = m0*m0    - m2*m2*hy;
    float k1 = m0*m0*m3 - m1*m2*hy + m0*ra*(rr*m2*1.0        );
    float k0 = m0*m0*m5 - m1*m1*hy + m0*ra*(rr*m1*2.0 - m0*ra);
    float h = k1*k1 - k2*k0;
    if( h<0.0 ) return vec4(MAX_DIST); //no intersection
    float t = (-k1-sqrt(h))/k2;
    float y = m1 + t*m2;
    if( y<0.0 || y>m0 ) return vec4(MAX_DIST); //no intersection
    return vec4(t, normalize(m0*(m0*(oa+t*rd)+rr*ba*ra)-ba*hy*y));
}

float itorus(vec3 pos,vec3 axe,float degree,  in vec2 tor ,out vec3 ro,out vec3 rd)
{
    mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;
    rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    float po = 1.0;
    float Ra2 = tor.x*tor.x;
    float ra2 = tor.y*tor.y;
    float m = dot(ro,ro);
    float n = dot(ro,rd);
    float k = (m + Ra2 - ra2)/2.0;
    float k3 = n;
    float k2 = n*n - Ra2*dot(rd.xy,rd.xy) + k;
    float k1 = n*k - Ra2*dot(rd.xy,ro.xy);
    float k0 = k*k - Ra2*dot(ro.xy,ro.xy);
    
    if( abs(k3*(k3*k3-k2)+k1) < 0.01 )
    {
        po = -1.0;
        float tmp=k1; k1=k3; k3=tmp;
        k0 = 1.0/k0;
        k1 = k1*k0;
        k2 = k2*k0;
        k3 = k3*k0;
    }
    
    float c2 = k2*2.0 - 3.0*k3*k3;
    float c1 = k3*(k3*k3-k2)+k1;
    float c0 = k3*(k3*(c2+2.0*k2)-8.0*k1)+4.0*k0;
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;
    float Q = c2*c2 + c0;
    float R = c2*c2*c2 - 3.0*c2*c0 + c1*c1;
    float h = R*R - Q*Q*Q;
    
    if( h>=0.0 )  
    {
        h = sqrt(h);
        float v = sign(R+h)*pow(abs(R+h),1.0/3.0); // cube root
        float u = sign(R-h)*pow(abs(R-h),1.0/3.0); // cube root
        vec2 s = vec2( (v+u)+4.0*c2, (v-u)*sqrt(3.0));
        float y = sqrt(0.5*(length(s)+s.x));
        float x = 0.5*s.y/y;
        float r = 2.0*c1/(x*x+y*y);
        float t1 =  x - r - k3; t1 = (po<0.0)?2.0/t1:t1;
        float t2 = -x - r - k3; t2 = (po<0.0)?2.0/t2:t2;
        float t = 1e20;
        if( t1>0.0 ) t=t1;
        if( t2>0.0 ) t=min(t,t2);
        return t;
    }
    
    float sQ = sqrt(Q);
    float w = sQ*cos( acos(-R/(sQ*Q)) / 3.0 );
    float d2 = -(w+c2); if( d2<0.0 ) return MAX_DIST;
    float d1 = sqrt(d2);
    float h1 = sqrt(w - 2.0*c2 + c1/d1);
    float h2 = sqrt(w - 2.0*c2 - c1/d1);
    float t1 = -d1 - h1 - k3; t1 = (po<0.0)?2.0/t1:t1;
    float t2 = -d1 + h1 - k3; t2 = (po<0.0)?2.0/t2:t2;
    float t3 =  d1 - h2 - k3; t3 = (po<0.0)?2.0/t3:t3;
    float t4 =  d1 + h2 - k3; t4 = (po<0.0)?2.0/t4:t4;
    float t = 1e20;
    if( t1>0.0 ) t=t1;
    if( t2>0.0 ) t=min(t,t2);
    if( t3>0.0 ) t=min(t,t3);
    if( t4>0.0 ) t=min(t,t4);
    return t;
}

float capIntersect(vec3 pos,vec3 axe,float degree,  in vec3 pa, in vec3 pb, in float ra ,out vec3 ro,out vec3 rd )
{
    mat4 translatem=translate(pos);
    mat4 rotatem=rotation(axe,degree);
    ro=(rotatem*translatem*INVMODEL*vec4(cameraPosition,1.0)).xyz;


     rd=normalize((rotatem* translatem*INVMODEL*vec4(vec3(uv,cameraPosition.z-1),1.0)).xyz-ro);
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);
    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - ra*ra*baba;
    float h = b*b - a*c;
    if( h >= 0.0 )
    {
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        // body
        if( y>0.0 && y<baba ) return t;
        // caps
        vec3 oc = (y <= 0.0) ? oa : ro - pb;
        b = dot(rd,oc);
        c = dot(oc,oc) - ra*ra;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return -1.0;
}

vec3 capNormal( in vec3 pos, in vec3 a, in vec3 b, in float r )
{
    vec3  ba = b - a;
    vec3  pa = pos - a;
    float h = clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0);
    return (pa - h*ba)/r;
}

vec3 nCylinder( in vec3 p, in vec3 a, in vec3 b, in float ra )
{
    vec3  pa = p - a;
    vec3  ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);
    return (pa - ba*paba/baba)/ra;
}

vec3 nTorus( in vec3 pos, vec2 tor )
{
	return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}

vec3 sphNormal( in vec3 pos, in vec3 centre )
{
    return normalize(pos-centre);
}
    vec3 ambientColor = vec3(0.1, 0.1, 0.1);
    vec3 lightPosition = vec3(0, 5, 5.0); 
    vec3 lightColor = vec3(1.0, 1.0, 1.0);    
    float specularStrength = 0.5; 



vec4 shading(in vec3 pos,in vec3 nor, in vec3 materialColor  ){


            

    vec3 ambient = ambientColor * materialColor;

    vec3 lightDir = normalize(lightPosition - pos);
    float diff = max(dot(nor, lightDir), 0.0);
    vec3 diffuse = diff * lightColor*materialColor;


    vec3 viewDir = normalize(cameraPosition - pos);
    vec3 reflectDir = reflect(-lightDir, nor);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = specularStrength * spec * lightColor;

  
    vec3 lighting = ambient + diffuse + specular;

 
  return vec4(lighting, 1.0);

}

void main()
{
  vec3 ro;
  vec3 rd;

  vec3 nor;
  vec3 pos;
 
 vec3 centre=vec3(4,0,0);
  vec2 res=isphere(centre,vec3(0,1,0),0,0.5,ro,rd); 
  float hit=min(res.x,res.y);
  if(hit>0&&hit<MAX_DIST)
  {
    pos = ro+hit*rd;
    nor = sphNormal(pos,centre);

    FragColor = shading(pos,nor,vec3(1, 0.0, 0.0));
  }  
    else    
    {
      vec3 a=vec3(0,0,0);
      vec3 b=vec3(0,1,0);
      float ra=0.1;
      vec4 res=icylinder(vec3(0.5,0,0),vec3(0,0,1),90,a,b,ra,ro,rd);
      hit=res.x;
        
      if(hit>0&&hit<MAX_DIST)
      {
        
        pos = ro+hit*rd;
        nor =res.yzw;  //nCylinder(pos, a,b,ra );
        FragColor =FragColor = shading(pos,nor,vec3(0.0, 1.0, 0.0));
      }
        
        else    
        {       
           res=coneIntersect(vec3(2,0,0),vec3(1,0,0),10,vec3(0,0,0),vec3(0,1 ,0),0.5,0.1,ro,rd);
           hit=res.x;         
           if(hit>0&&hit<MAX_DIST)
            {
                pos = ro+hit*rd;
                nor =  res.yzw;      
                FragColor = FragColor = shading(pos,nor,vec3(0.0, 0.0, 1.0));  
            }               
            else    
              {
                vec2 Rr=vec2(1.0,0.5);
                float hit=itorus(vec3(-4,0,0),vec3(0,1,0),20,Rr,ro,rd);        
                if(hit>0&&hit<MAX_DIST)
                {
                 pos = ro+hit*rd;
                 nor= nTorus(pos, Rr );
                 FragColor = vec4(1.0, 1.0, 0.0, 1.0);  
                }
                      
                  else    
                  {
                    vec3 a=vec3(0,0,0);
                    vec3 b=vec3(1,0,0);
                    float ra=0.5;
                    float hit=capIntersect(vec3(-2,0,0),vec3(0,0,1),0,a,b,ra,ro,rd );        
                    if(hit>0&&hit<MAX_DIST)
                    {
                        pos = ro+hit*rd;
                        nor = capNormal(pos, a,b,ra );
                        FragColor = vec4(0.0, 1.0, 0.0, 1.0);  
                    }
                    
                      else
                      {
                       nor=vec3(0,1,0);
                       hit=iplane(vec3(0,0,0),vec3(1,0,0),0,vec3(0,0,0),nor,ro,rd); 
                       if(hit>0&&hit<MAX_DIST)
                       {
                          
                           pos=ro+hit*rd;                        
                            FragColor = shading(pos,nor,vec3(0.5, 0.5, 0.5));
                       }
                                         
                       else
                        FragColor = vec4(0.0, 0.0, 0, 1.0);  
                        }

                  }
              }


        }
         
      
    }


} 

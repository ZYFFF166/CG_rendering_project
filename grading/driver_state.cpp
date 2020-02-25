#include "driver_state.h"
#include <cstring>
#include <vector>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    unsigned int area = width * height;
    state.image_color = new pixel[area];
    state.image_depth = new float[area];
    for (unsigned int i = 0; i < area; ++i)
    {
        state.image_color[i] = make_pixel(0, 0, 0);
    }
    for (size_t i = 0; i < area; ++i)
    {
        state.image_depth[i] = -200.0;
    }   

}

	//std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    switch (type)
    {
        case render_type::triangle:
        {
            data_geometry arr[3];
            const data_geometry * ptr[3] = { &arr[0], &arr[1], &arr[2] };
            for (int vert = 0; vert < state.num_vertices * state.floats_per_vertex; vert += 3 * state.floats_per_vertex)
            {
                for (int v = 0; v < 3; ++v)
                {
                    arr[v].data = state.vertex_data + vert + v * state.floats_per_vertex;
                }
                //rasterize_triangle(state, ptr);
                clip_triangle(state, ptr, 0);
            }
            break;
        }
        case render_type::indexed:
        {
            data_geometry center[3];
            const data_geometry * ptr[3] = { &center[0],&center[1],&center[2] };
            size_t max = 3 * state.num_triangles;
            
            for (size_t i = 0; i < max; i += 3)
            {
                center[0].data = state.vertex_data + *(state.index_data + i) * 3;
                center[1].data = state.vertex_data + *(state.index_data + i + 1) * 3;
                center[2].data = state.vertex_data + *(state.index_data + i + 2) * 3;
                //rasterize_triangle(state, ptr);
                clip_triangle(state, ptr, 0);
            }
            break;
        }
        case render_type::fan:
        {
            data_geometry center[3];
            const data_geometry * ptr[3] = { &center[0],&center[1],&center[2] };
            size_t max = state.num_vertices * state.floats_per_vertex;
            size_t add = 3 * state.floats_per_vertex;
            
            center[0].data = state.vertex_data;
            center[1].data = state.vertex_data + state.floats_per_vertex;
            center[2].data = state.vertex_data + 2 * state.floats_per_vertex;
            //rasterize_triangle(state, ptr);
            clip_triangle(state, ptr, 0);
            for (size_t i = add; i < max; i += state.floats_per_vertex)
            {
                center[1].data = center[2].data;
                center[2].data = state.vertex_data + i;
                //rasterize_triangle(state, ptr);
                clip_triangle(state, ptr, 0);
            }
            break;
        }
        case render_type::strip:
        {
            data_geometry center[3];
            const data_geometry * ptr[3] = { &center[0],&center[1],&center[2] };
            size_t max = state.num_vertices * state.floats_per_vertex;
            size_t add = 3 * state.floats_per_vertex;
            
            center[0].data = state.vertex_data;
            center[1].data = state.vertex_data + state.floats_per_vertex;
            center[2].data = state.vertex_data + 2 * state.floats_per_vertex;
            //rasterize_triangle(state, ptr);
            clip_triangle(state, ptr, 0);
            
            for (size_t i = add; i < max; i += state.floats_per_vertex)
            {
                center[0].data = center[1].data;
                center[1].data = center[2].data;
                center[2].data = state.vertex_data + i;
                //rasterize_triangle(state, ptr);
                clip_triangle(state, ptr, 0);
            }
            break;
        }
        default: {}
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    data_vertex v;
    data_geometry out[3];
    const data_geometry * p_out[3] = {&out[0] ,&out[1], &out[2]};
    
    for(unsigned int i = 0; i < 3; ++i)
    {
        v.data = in[i]->data;
        out[i].data = v.data;
        out[i].gl_Position = in[i]->gl_Position;
        if (face == 0)
        {
            state.vertex_shader(v, out[i], state.uniform_data);
        }
    }
    vec3 a = {out[0].gl_Position[0] , out[0].gl_Position[1] , out[0].gl_Position[2]};
    vec3 b = {out[1].gl_Position[0] , out[1].gl_Position[1] , out[1].gl_Position[2]};
    vec3 c = {out[2].gl_Position[0] , out[2].gl_Position[1] , out[2].gl_Position[2]};

    vec3 point;
    vec3 norm;
    
    std::vector<vec3> clipped_in, clipped_out;
    std::vector<data_geometry> holder;
    std::vector<data_geometry> goner;
    vec3 clipped_v1, clipped_v2, u1,v1,u2,v2;
	if(face==0)
    {
        point ={0,0,1};
        norm = {0,0,-1};   
    }
     else if(face==1)
    {
        point ={0,0,-1};
        norm = {0,0,1};   
    }
 	else if(face==2)
    {
        point ={-1,0,0};
        norm = {1,0,0};
    }
      else if(face==3)
    {
        point ={1,0,0};
        norm = {-1,0,0};
    }
     else if(face==4)
    {
        point ={0,-1,0};
        norm = {0,1,0}; 
    } 
    else if(face==5)
    {
        point ={0,1,0};
        norm = {0,-1,0};  
    }
    else if(face==6)
    {
        rasterize_triangle(state, p_out);
        return;
    }
    
    
    if(dot(norm,(point-a)) <= out[0].gl_Position[3])
    {
        clipped_in.push_back(a);
        holder.push_back(out[0]);
    }
    else 
    {
        clipped_out.push_back(a);
        goner.push_back(out[0]);
    }
    
    if(dot(norm,(point-b)) <= out[1].gl_Position[3])
    {
        clipped_in.push_back(b);
        holder.push_back(out[1]);
    }
    else 
    {
        clipped_out.push_back(b);
        goner.push_back(out[1]);
    }
    
    if(dot(norm,(point-c)) <= out[2].gl_Position[3])
    {
        clipped_in.push_back(c);
        holder.push_back(out[2]);
    }
    else 
    {
        clipped_out.push_back(c);
        goner.push_back(out[2]);
    }
    
    
    data_geometry out_1[3];
    const data_geometry * p_out_1[3] = {&out_1[0] ,&out_1[1] ,&out_1[2]} ;


    if(clipped_in.size() == 0)
    { 
    	return;
    }
    if(clipped_in.size() == 1)
    {
  		//out_1[1]
        out_1[0].data = holder[0].data;
        out_1[0].gl_Position = holder[0].gl_Position;

        u1 = clipped_out[0] - clipped_in[0];
        v1 = clipped_in[0] - point;

        clipped_v1 = clipped_in[0] + dot(norm, v1) / dot(norm,u1) * u1;
        
        out_1[1].data = goner[0].data;
        out_1[1].gl_Position = goner[0].gl_Position;       
        out_1[1].gl_Position[0] = clipped_v1[0];
        out_1[1].gl_Position[1] = clipped_v1[1];
        out_1[1].gl_Position[2] = clipped_v1[2];
        
		//out_1[2]
        u2 = clipped_out[1] - clipped_in[0];
        v2 = clipped_in[0] - point;
        
        clipped_v2 = clipped_in[0] + dot(norm, v2) / dot(norm,u2) * u2;
        
        out_1[2].data = goner[1].data;
        out_1[2].gl_Position = goner[1].gl_Position;
        
        out_1[2].gl_Position[0] = clipped_v2[0];
        out_1[2].gl_Position[1] = clipped_v2[1];
        out_1[2].gl_Position[2] = clipped_v2[2];
        
        clip_triangle(state, p_out_1, face+1);
    }
    if(clipped_in.size() == 2)
    {
        
    }
    
    clip_triangle(state, p_out, face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
	int walk = 0;
    float alpha, beta, gamma;  
    size_t w = state.image_width;
    size_t h = state.image_height;
    vec2 p;
    vec3 p_3d;
    float a_x,a_y,b_x,b_y,c_x,c_y;
    a_x = (in[0]->gl_Position[0] /in[0]->gl_Position[3] + 1) * 0.5 * w;
    a_y = (in[0]->gl_Position[1] /in[0]->gl_Position[3] + 1) * 0.5 * h;
    b_x = (in[1]->gl_Position[0] /in[1]->gl_Position[3] + 1) * 0.5 * w; 
    b_y = (in[1]->gl_Position[1] /in[1]->gl_Position[3] + 1) * 0.5 * h;
    c_x = (in[2]->gl_Position[0] /in[2]->gl_Position[3] + 1) * 0.5 * w;
    c_y = (in[2]->gl_Position[1] /in[2]->gl_Position[3] + 1) * 0.5 * h;
    vec2 a = vec2(a_x,a_y);  
    vec2 b = vec2(b_x,b_y); 
    vec2 c = vec2(c_x,c_y);
    
    for (size_t i = 0; i < h; ++i)
    {
        for (size_t j = 0; j < w; ++j)
        {
            p = {float(j + 0.5), float(i + 0.5)};
			float area_abc = 0.5 * ((b[0]*c[1]-c[0]*b[1])-(a[0]*c[1]-c[0]*a[1])+(a[0]*b[1]-b[0]*a[1]));
			float area_pbc = 0.5 * ((b[0]*c[1]-c[0]*b[1])-(p[0]*c[1]-c[0]*p[1])+(p[0]*b[1]-b[0]*p[1]));  
			
			float area_apc = 0.5 * ((p[0]*c[1]-c[0]*p[1])-(a[0]*c[1]-c[0]*a[1])+(a[0]*p[1]-p[0]*a[1])); 
			float area_abp = 0.5 * ((b[0]*p[1]-p[0]*b[1])-(a[0]*p[1]-p[0]*a[1])+(a[0]*b[1]-b[0]*a[1]));
			alpha = area_pbc/area_abc;
			beta = area_apc/area_abc;
			gamma = area_abp/area_abc;
            
            p_3d[0] = p[0];
            p_3d[1] = p[1];
            p_3d[2] = -(in[0]->gl_Position[2]/in[0]->gl_Position[3]*alpha + in[1]->gl_Position[2]/in[1]->gl_Position[3]*beta + in[2]->gl_Position[2]/in[2]->gl_Position[3]*gamma);
            
            if (alpha >= 0 && beta >= 0 && gamma >= 0)
            {
               	float d_out[state.floats_per_vertex];
                data_output output;
                //data_fragment frag;
                for (int k = 0; k < state.floats_per_vertex; ++k)
                {
                    if(state.interp_rules[k] == interp_type::flat)
                    {
                        d_out[k]= in[0]->data[k];
                    }
                    else if(state.interp_rules[k] == interp_type::smooth)
                    {
                        d_out[k]= ((alpha * in[0]->data[k]/in[0]->gl_Position[3]) + (beta * in[1]->data[k]/in[1]->gl_Position[3]) + (gamma * in[2]->data[k]/in[2]->gl_Position[3]))
                        /((alpha/in[0]->gl_Position[3]) + (beta/in[1]->gl_Position[3]) + (gamma/in[2]->gl_Position[3]));
                    }
                    else if(state.interp_rules[k] == interp_type::noperspective)
                    {
                        d_out[k]= (alpha * in[0]->data[k]) + (beta * in[1]->data[k]) + (gamma * in[2]->data[k]);
                    }
                    
                }
                
                if(p_3d[2] > state.image_depth[walk] && p_3d[2] <= 1  && p_3d[2] >= -1)
                {
                    state.image_depth[walk] = p_3d[2];
                    state.fragment_shader({d_out}, output, state.uniform_data);
                    state.image_color[walk] = make_pixel(255*output.output_color[0],255*output.output_color[1],255*output.output_color[2]);
                }
                    
            }
            ++walk;
        }
    }
}














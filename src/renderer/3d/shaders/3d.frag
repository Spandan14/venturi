#version 410 core

in vec2 uv;
out vec4 frag_color;

uniform sampler3D density_tex;

uniform vec3 camera_pos;
uniform mat4 inv_view;
uniform mat4 inv_proj;

uniform vec3 volume_min;
uniform vec3 volume_max;

uniform float step_size;
uniform float absorption;
uniform int max_steps;

vec3 compute_ray_direction() {
    vec2 ndc = uv * 2.0 - 1.0;

    vec4 clip = vec4(ndc, -1.0, 1.0);
    vec4 view = inv_proj * clip;
    view /= view.w;

    vec3 world_dir = normalize((inv_view * vec4(view.xyz, 0.0)).xyz);

    return world_dir;
}

bool intersect_aabb(vec3 ray_origin, vec3 ray_dir, vec3 box_min, vec3 box_max, out float t_near, out float t_far) {
    vec3 inv_d = 1.0 / ray_dir;
    vec3 t_min = (box_min - ray_origin) * inv_d;
    vec3 t_max = (box_max - ray_origin) * inv_d;

    vec3 t1 = min(t_min, t_max);
    vec3 t2 = max(t_min, t_max);

    t_near = max(max(t1.x, t1.y), t1.z);
    t_far = min(min(t2.x, t2.y), t2.z);

    return t_far >= max(t_near, 0.0);
}

float sample_density(vec3 pos) {
    vec3
    tex_coord_3d = (pos - volume_min) / (volume_max - volume_min);
    return texture(density_tex, tex_coord_3d).r;
}

void main() {
    vec3 ray_orig = camera_pos;
    vec3 ray_dir = compute_ray_direction();

    float t_entry, t_exit;
    if (!intersect_aabb(ray_orig, ray_dir, volume_min, volume_max, t_entry, t_exit)) {
        frag_color = vec4(0.0);
        return;
    }

    t_entry = max(t_entry, 0.0);
    vec3 color = vec3(0.0);
    float alpha = 0.0;

    float t = t_entry;
    for (int i = 0; i < max_steps; ++i) {
        if (t > t_exit || alpha >= 0.999) {
            break;
        }

        vec3 pos = ray_orig + t * ray_dir;
        float density = sample_density(pos);

        float a = 1.0 - exp(-density * absorption);
        vec3 c = vec3(density);

        color += (1.0 - alpha) * a * c;
        alpha += (1.0 - alpha) * a;

        t += step_size;
    }

    frag_color = vec4(color, alpha);
}

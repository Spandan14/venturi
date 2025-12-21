#version 410 core

in vec2 uv;
out vec4 frag_color;

uniform sampler3D density_tex;

uniform vec3 camera_pos;
uniform mat4 inv_view;
uniform mat4 inv_proj;

uniform vec3 volume_min;
uniform vec3 volume_max;
uniform ivec3 grid_size;

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

        vec3 grid_pos = (pos - volume_min) / (volume_max - volume_min) * vec3(grid_size);
        vec3 frac = fract(grid_pos);
        vec3 dist_to_face = min(frac, 1.0 - frac);
        float border_width = 0.00;
        vec3 border_color = vec3(0.0, 0.0, 0.0);

        int lo_borders = 0;
        int hi_borders = 0;
        if (pos.x < volume_min.x + 0.5) {
            lo_borders += 1;
        }
        if (pos.x > volume_max.x - 0.5) {
            hi_borders += 1;
        }
        if (pos.y < volume_min.y + 0.5) {
            lo_borders += 1;
        }
        if (pos.y > volume_max.y - 0.5) {
            hi_borders += 1;
        }
        if (pos.z < volume_min.z + 0.5) {
            lo_borders += 1;
        }
        if (pos.z > volume_max.z - 0.5) {
            hi_borders += 1;
        }

        if (lo_borders >= 2) {
            border_width = 0.05;
            border_color = vec3(0.0, 1.0, 0.0);
        }

        if (hi_borders >= 2) {
            border_width = 0.05;
            border_color = vec3(1.0, 0.0, 0.0);
        }

        float border = step(dist_to_face.x, border_width) + step(dist_to_face.y, border_width) + step(dist_to_face.z, border_width);

        border = clamp(border, 0.0, 1.0);

        float density = sample_density(pos);

        // guard
        if (density < 1e-3) {
            t += step_size;
            continue;
        }

        float a = 1.0 - exp(-density * absorption);
        vec3 c = mix(vec3(density), border_color, border);
        // c = vec3(density);

        color += (1.0 - alpha) * a * c;
        alpha += (1.0 - alpha) * a;

        t += step_size;

        // color = grid_pos;
        // alpha = 1.0;
    }

    frag_color = vec4(color, alpha);
    // frag_color = vec4(t, t_entry, t_entry, 1.0);
    // frag_color = vec4(ray_dir.x, ray_dir.y, ray_dir.z, 1.0);
}

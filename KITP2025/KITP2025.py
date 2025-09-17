from manim import *
from manim_slides import Slide
from numpy import cos, sin, pi
from MF_Tools import TransformByGlyphMap

class KITP2025(Slide):
    def construct(self):
        self.camera.background_color = LIGHTER_GRAY
        title = VGroup(
            Paragraph("""Strongly Correlated Physics in Rhombohedral
            Few Layer Graphene""", font="Times New Roman", line_spacing=1, color=BLACK, alignment="center"),
            Paragraph("Ian Sackin, PI: Andrea Young, Mentor: Ludwig Holleis", color=LOGO_BLUE).scale(0.4)
        ).arrange(DOWN, buff=.7)
        title.add(SurroundingRectangle(title[0], color=LOGO_BLACK, buff=0.3))

        self.play(DrawBorderThenFill(title))
        self.wait()

        self.next_slide()

        self.play(
            title[0][0][:25].animate.scale(0.6).move_to([-3.5,2.7,0]).set(weight=BOLD),
            FadeOut(title[0][0][25:]), FadeOut(title[0][1]),
            FadeOut(VGroup(title[1:]))
        )
        self.wait()

        vague_def = Tex("electron-electron interaction $>$ KE").scale(1.2)
        effects = Tex("chiral superconductivity$^{[1]}$, FQH effect$^{[2]}$, etc.")
        coulomb_energy = MathTex(r"U \sim \frac{e^2}{\epsilon} \frac{1}{r} ").scale(1.2)
        energy_scale = Tex("increase electron density $n_e$").scale(1.2)
        VGroup(vague_def, effects, coulomb_energy, energy_scale).arrange(DOWN, buff=0.5).set_color(BLACK)

        self.next_slide()
        self.play(Write(vague_def))
        self.wait()

        self.next_slide()
        self.play(Write(effects))
        self.wait()

        self.next_slide()
        self.play(Write(coulomb_energy))
        self.wait()

        self.next_slide()
        self.play(Indicate(coulomb_energy[0][8]), run_time=3, color=RED)
        self.wait()

        self.next_slide()
        self.play(Write(energy_scale))
        self.wait()

        self.next_slide()
        self.play(
            Unwrite(VGroup(vague_def, coulomb_energy, effects)),
            energy_scale.animate.scale(0.6).move_to([-3.7,0,0])
        )

        self.wait()

        

        ## Band Physics
        self.next_slide()
        band_physics = Tex("Band Physics",color=BLACK).move_to([-4,2.7,0])
        self.play(
            title[0][0][:25].animate.shift(UP*3),
            Write(band_physics)
        )
        self.wait()

        ax = Axes(
            axis_config={"include_ticks": False, "include_numbers": False, "include_tip":False},
            x_length=8
        ).set_color(BLACK).shift(RIGHT*1)
        ax.x_axis.shift(DOWN*2)
        labels = VGroup(
            MathTex(r"\epsilon", color=BLACK).next_to(ax.y_axis.get_top(), RIGHT),
            MathTex("k", color=BLACK).next_to(ax.x_axis.get_right(), UP)
        )

        self.play(Write(ax), Write(labels))
        self.wait()
        
        self.next_slide()

        cband = ax.plot(lambda x: (x/2)**2+.5, x_range=[-PI,PI]).set_color(BLUE)
        vband = ax.plot(lambda x: -(x/3)**2 - 1.5, x_range=[-PI,PI]).set_color(BLUE)

        self.play(Write(cband), Write(vband))
        self.wait()

        # Discrete k-points for electrons
        k_vals = np.linspace(-PI, PI, 21)

        # Fermi level tracker
        fermi_level = ValueTracker(0.0)

        # Fermi level dashed line
        fermi_line = always_redraw(
            lambda: ax.plot( 
                lambda x: fermi_level.get_value(),
                x_range=[-PI, PI]
            ).set_color(RED)
        )
        fermi_line_label = always_redraw(lambda: MathTex(r"\epsilon_F").next_to(fermi_line, RIGHT).set_color(BLACK))

        # Electrons as dots, update depending on Fermi level

        dots = always_redraw(
            lambda: VGroup(*[
                Dot(ax.coords_to_point(k, vband.underlying_function(k)), radius=0.05, color=BLACK)
                for k in k_vals if vband.underlying_function(k) <= fermi_level.get_value()
            ] + [
                Dot(ax.coords_to_point(k, cband.underlying_function(k)), radius=0.05, color=BLACK)
                for k in k_vals if cband.underlying_function(k) <= fermi_level.get_value()
            ])
        )

        dotsv = always_redraw(
            lambda: VGroup(*[Dot(ax.coords_to_point(k, vband.underlying_function(k)), radius=0.05, color=DARK_BLUE)
                for k in k_vals if vband.underlying_function(k) <= fermi_level.get_value()])
        )
        dotsc = always_redraw( 
                lambda: VGroup(*[Dot(ax.coords_to_point(k, cband.underlying_function(k)), radius=0.05, color=DARK_BLUE)
                for k in k_vals if cband.underlying_function(k) <= fermi_level.get_value()])
        )


        self.next_slide()
        self.play(Create(fermi_line), Create(fermi_line_label), FadeIn(dotsv, dotsc))

        # Animate moving the Fermi level down into the valence band
        self.next_slide()
        self.play(fermi_level.animate.set_value(-2.0), run_time=3)
        # Animate moving the Fermi level up into the conduction band
        self.next_slide()
        self.play(fermi_level.animate.set_value(2.0), run_time=3)
        self.wait()


        flat_band = ax.plot(lambda x: (x/3)**6+.5, x_range=[-PI,PI]).set_color(BLUE)
        dotsf = always_redraw( 
                lambda: VGroup(*[Dot(ax.coords_to_point(k, flat_band.underlying_function(k)), radius=0.05, color=DARK_BLUE)
                for k in k_vals if flat_band.underlying_function(k) <= fermi_level.get_value()])
        )

        self.next_slide()
        self.play(fermi_level.animate.set_value(0.5))
        self.play(
            ReplacementTransform(cband, flat_band),
        )
        self.remove(dotsc)
        self.add(dotsf.update())
        self.wait()

        self.next_slide()
        self.play(fermi_level.animate.set_value(0.6), run_time=6, rate_func=rush_from)
        self.wait()

        self.next_slide()
        self.play(*[FadeOut(mob) for mob in self.mobjects])

        ## Graphene
        # --- Parameters ---
        hex_r = 1                      # circumradius of hexagon
        hex_corners = [
            [hex_r * np.cos(theta), hex_r * np.sin(theta), 0]
            for theta in np.linspace(0, 2 * np.pi, 7)  # 7 points closes the polygon
        ]
        true_r = hex_r * np.sqrt(3) / 2  # apothem

        hexagon = Polygon(*hex_corners)
        hex_connections = VGroup()

        dx = np.sqrt(3) * true_r
        dy = 2 * true_r  

        # --- Tile hexagons ---
        for i in range(-5, 6):
            for j in range(-5, 6):
                x = j * dx
                y = i * dy
                if j % 2 == 1:
                    y += true_r
                hex_connections.add(
                    hexagon.copy().move_to([x, y, 0])
                )

        hex_connections.set_stroke(opacity=0.3)

        # --- Collect unique lattice points ---
        lattice_points = []
        eps = 1e-6
        def unique(pt):
            for q in lattice_points:
                if np.linalg.norm(pt - q) < eps:
                    return False
            return True

        for hex in hex_connections:
            c = hex.get_center()
            for angle in range(6):
                pt = c + RIGHT * hex_r
                pt = rotate_vector(pt - c, angle * PI/3) + c
                if unique(pt):
                    lattice_points.append(pt)

        # --- Build atoms ---
        hex_lattice = VGroup()
        for k, pt in enumerate(lattice_points):
            color = BLUE if k % 2 == 0 else BLUE  # simple sublattice coloring
            hex_lattice.add(Dot(pt, radius=0.07, color=color))

        # --- Center everything around origin ---
        # pick the dot closest to origin
        center_point = min(lattice_points, key=lambda p: np.linalg.norm(p))
        shift_vec = -center_point
        hex_connections.shift(shift_vec)
        hex_lattice.shift(shift_vec)

        A_lattice = VGroup(hex_connections, hex_lattice)
        B_lattice = A_lattice.copy().set_color(ORANGE).shift(hex_r * LEFT)
        C_lattice = B_lattice.copy().set_color(PURPLE).shift(hex_r * np.array([cos(pi/3), sin(pi/3), 0]))

        axes = Axes(
            x_range=[-2,2],
            y_range=[-2,2],
            axis_config={"include_ticks": False, "include_numbers": False, "include_tip":False},
            x_length=4, y_length=6
        ).to_edge(RIGHT, buff=1.5).set_z_index(101).set_color(BLACK)
        axes_labels = VGroup(
            MathTex(r"\epsilon", color=BLACK).next_to(axes.y_axis.get_top(), RIGHT),
            MathTex("k", color=BLACK).next_to(axes.x_axis.get_right(), UP)
        ).set_z_index(101)
        axes_background = BackgroundRectangle(axes, buff=0.2).set_z_index(90).set_color(LOGO_WHITE).set_opacity(1)

        # MLG Dirac cone: E ~ ±k
        mlg_band = VGroup(
            axes.plot(lambda x: np.abs(x), color=BLUE),
            axes.plot(lambda x: -np.abs(x), color=BLUE)
        ).set_z_index(91)

        # Bernal bilayer: E ~ ±k^2
        bilayer_band = VGroup(
            axes.plot(lambda x: x**2/2, color=ORANGE),
            axes.plot(lambda x: -x**2/2, color=ORANGE)
        ).set_z_index(92)

        # ABA multilayer: bilayer + monolayer
        aba_band = VGroup(
            axes.plot(lambda x: x, color=ORANGE),
            axes.plot(lambda x: -x, color=ORANGE),
            axes.plot(lambda x: x**2/2, color=ORANGE),
            axes.plot(lambda x: -x**2/2, color=ORANGE)
        ).set_z_index(93)

        # ABC: cubic dispersion E ~ ±k^3
        abc3_band = VGroup(
            axes.plot(lambda x: np.abs(x**3)/2, x_range=[-1.6,1.6], color=GREEN),
            axes.plot(lambda x: -np.abs(x**3)/2, x_range = [-1.6,1.6], color=GREEN)
        ).set_z_index(94)

        # D-field Splitting E ~ ±k^2N
        split_band = VGroup(
            axes.plot(lambda x: x**6/2 + 0.5, x_range=[-1.2,1.2], color=GREEN),
            axes.plot(lambda x: -x**6/2 - 0.5, x_range = [-1.2,1.2], color=GREEN)
        ).set_z_index(94)

        # k^N generalization for ABC
        def abc_band(n, color=RED):
            return VGroup(
                axes.plot(lambda x: x**n/2, color=color, x_range=[-1.5,1.5]),
                axes.plot(lambda x: -x**n/2, color=color, x_range=[-1.5,1.5])
            ).set_z_index(95)



        # --- Add to scene ---
        self.next_slide()
        mlg_title = Tex("Mono-Layer Graphene").move_to([-4,2.7,0]).set_color(BLACK)
        mlg_title = VGroup( BackgroundRectangle(mlg_title, buff=0.2).set_color(GRAY).set_opacity(0.7), mlg_title).set_z_index(100)
        self.play(DrawBorderThenFill(mlg_title))
        self.wait()

        self.play(DrawBorderThenFill(A_lattice), run_time =2)
        self.wait()

        ## MLG Band Structure
        self.next_slide()
        self.play(Write(axes), Write(axes_labels), FadeIn(axes_background))
        self.wait()

        self.next_slide()
        self.play(Create(mlg_band))
        self.wait()


        self.next_slide()
        self.play(Unwrite(mlg_title[1]))
        self.wait()

        ## Bernal Stacking

        bernal_title = Tex("ABA/Bernal Stacking").move_to([-4,2.7,0]).set_color(BLACK)
        bernal_title = VGroup( BackgroundRectangle(bernal_title, buff=0.2).set_color(GRAY).set_opacity(0.7), bernal_title).set_z_index(100)
        stable = Tex("Most Stable").scale(0.8).set_color(BLACK).next_to(bernal_title, DOWN)
        stable = VGroup( BackgroundRectangle(stable, buff=0.2).set_color(GRAY).set_opacity(0.7), stable).set_z_index(100)

        self.next_slide()
        self.play(DrawBorderThenFill(B_lattice))
        self.wait()

        self.next_slide()
        self.play(
            ReplacementTransform(mlg_title[0], bernal_title[0]),
            Write(bernal_title[1]), Write(stable)
        )
        self.wait()

        self.next_slide()
        self.play(
            ReplacementTransform(mlg_band[0], bilayer_band[0]),
            ReplacementTransform(mlg_band[1], bilayer_band[1])
        )
        self.wait()

        aba_cite = Paragraph("""[4]: M. Koshino and E. McCann, Electronic properties of 
                             multilayer graphene, Phys. Rev. B 80, 165409 (2009)""").set_color(RED).scale(0.3).next_to(axes_background,DOWN).set_z_index(50)

        self.next_slide()
        self.play(FadeIn(aba_band), FadeOut(bilayer_band), FadeIn(aba_cite))
        self.wait()

        self.next_slide()
        self.play(DrawBorderThenFill(C_lattice))
        self.wait()

        rhombohedral_title = Tex("ABC/Rhombohedral Stacking").move_to([-3.5,2.7,0]).set_color(BLACK)
        rhombohedral_title = VGroup( BackgroundRectangle(rhombohedral_title, buff=0.2).set_color(GRAY).set_opacity(0.7), rhombohedral_title).set_z_index(100)
        metastable = Tex("Meta-Stable").scale(0.8).set_color(BLACK).next_to(rhombohedral_title, DOWN)
        metastable = VGroup( BackgroundRectangle(metastable, buff=0.2).set_color(GRAY).set_opacity(0.7), metastable).set_z_index(100)

        self.next_slide()
        self.play(
            ReplacementTransform(bernal_title[0], rhombohedral_title[0]),
            ReplacementTransform(bernal_title[1], rhombohedral_title[1]),
            ReplacementTransform(stable[0], metastable[0]),
            ReplacementTransform(stable[1], metastable[1])
        )
        self.wait()

        self.next_slide()
        self.play(ReplacementTransform(aba_band, abc3_band), FadeOut(aba_cite))
        self.wait()

        D0_dispersion = MathTex(r"\epsilon(k) \sim \pm |k|^N").scale(1.5).move_to([-3.2,-1.5,0]).set_color(BLACK).set_z_index(90)
        D0_dispersion_background = BackgroundRectangle(D0_dispersion, buff=0.2).set_color(GRAY).set_opacity(0.7)


        dispersion_cite = Paragraph("""[3]: Fan Zhang, Bhagawan Sahu, Hongki Min, and A. H. MacDonald, 
                                    Band structure of ABC-stacked graphene trilayers, Phys. Rev. B 82, 035409 (2010)"""
                              ).scale(0.3).set_color(RED).next_to(axes_background,DOWN).set_z_index(50)

        self.next_slide()
        self.play(FadeIn(D0_dispersion_background), Write(D0_dispersion), FadeIn(dispersion_cite))
        self.wait()

        D0_dispersion1 = MathTex(r"\epsilon(k) \sim \pm \sqrt{k^{2N}}").scale(1.5).move_to([-3.2,-1.5,0]).set_color(BLACK).set_z_index(90)
        D0_dispersion_background1 = BackgroundRectangle(D0_dispersion1, buff=0.2).set_color(GRAY).set_opacity(0.7)

        self.next_slide()
        self.play(
            TransformByGlyphMap(
                D0_dispersion, D0_dispersion1,
                ([7], [8]), ([9], [10]), ([6,8], FadeOut), (Write,[6,7,9])
            ),
            ReplacementTransform(D0_dispersion_background, D0_dispersion_background1)
        )
        self.wait()

        D_arrow = VGroup(
            Circle(radius=0.4).set_stroke(RED).set_fill(opacity=1,color=BLACK), Dot(color=RED)
        )
        D_arrow.add(MathTex(r"\vec{D}").set_color(BLACK).next_to(D_arrow[0], UL, buff=0)).scale(1.5).shift(UP*0.5)

        self.next_slide()
        self.play(Write(D_arrow))
        self.wait()

        Df_dispersion = MathTex(r"\epsilon(k) \sim \pm \sqrt{k^{2N} + \Delta^2}").scale(1.5).move_to([-3.2,-1.5,0]).set_color(BLACK).set_z_index(90)
        Df_dispersion_background = BackgroundRectangle(Df_dispersion, buff=0.2).set_color(GRAY).set_opacity(0.7)

        self.next_slide()
        self.play(
            TransformByGlyphMap(
                D0_dispersion1, Df_dispersion,
                (Write,[11,12,13])
            ),
            ReplacementTransform(D0_dispersion_background1, Df_dispersion_background)
        )
        self.wait()

        self.next_slide()
        self.play(
            ReplacementTransform(abc3_band[0], split_band[0]),
            ReplacementTransform(abc3_band[1], split_band[1])
        )
        self.wait()

        Df_dispersion1 = MathTex(r"\epsilon(k) \sim \pm k^{2N} \pm \Delta ").scale(1.5).move_to([-3.2,-1.5,0]).set_color(BLACK).set_z_index(90)
        Df_dispersion_background1 = BackgroundRectangle(Df_dispersion, buff=0.2).set_color(GRAY).set_opacity(0.7)

        self.play(
            ReplacementTransform(Df_dispersion_background, Df_dispersion_background1),
            TransformByGlyphMap(
                Df_dispersion, Df_dispersion1,
                ([6,7,13], FadeOut), ([12], [10]),
                ([8,9,10],[6,7,8]), ([11],[9])
            )
        )
        self.wait()

        self.next_slide()
        self.play(
            Unwrite(VGroup(A_lattice, B_lattice)), Unwrite(C_lattice, reverse=True),
            Unwrite(D_arrow), Unwrite(Df_dispersion1), Unwrite(Df_dispersion_background), Unwrite(dispersion_cite)
        )
        self.play(*[
            FadeOut(mob) for mob in self.mobjects
        ])

        self.wait()

        self.camera.background_color = LIGHTER_GRAY
        # --- Colors ---
        polymer_color = GREY
        bn_color = BLUE
        graphene_color = PINK
        gate_color = PURPLE
        substrate_color = PURPLE_A

        # --- Parameters ---
        x_left, x_right = -3, 3
        line_thickness = 16
        layer_spacing = 0.25
        pick_up_height = 2

        # --- Substrate ---
        substrate = Line([-5,0,0],[5,0,0], stroke_width=line_thickness).set_color(substrate_color)
        substrate.shift(DOWN*0.5)
        substrate_label = Tex("SiO$_2$ Substrate", font_size=24).next_to(substrate, RIGHT, buff=.2).set_color(BLACK)
        self.add(substrate)
        self.add(substrate_label)

        # Track stack top
        stack_top_y = 0.1

        # --- Helper function for layers + labels ---
        def create_layer(color, width_scale=1.0, label_text=""):
            layer = Line([x_left*width_scale,0,0],[x_right*width_scale,0,0], stroke_width=line_thickness).set_color(color)
            label = Text(label_text, font_size=24).next_to(layer, RIGHT, buff=0.2)
            return layer, label

        # --- Step 0: Polymer above ---
        polymer, polymer_label = create_layer(polymer_color, 1.7, "Polymer")
        polymer.shift(UP*3)
        polymer_label.shift(UP*3).set_color(BLACK)
        self.add(polymer, polymer_label)
        held_layers = [polymer, polymer_label]

        self.wait(0.5)

        flakes_info = [
            (bn_color, "hBN", 1.3),
            (graphene_color, "RG", 1),
            (bn_color, "hBN", 1.3),
            (gate_color, "BG", .7)
        ]

        for i, (color, label_text, width) in enumerate(flakes_info):

            self.next_slide()
            flake, flake_label = create_layer(color, width, label_text)
            flake_label.set_color(BLACK)
            self.play(FadeIn(VGroup(flake, flake_label)))

            polymer_target_y = stack_top_y + layer_spacing
            self.play(Group(*held_layers).animate.shift(DOWN*(polymer.get_bottom()[1]-polymer_target_y)), run_time=1)

            self.play(Group(*held_layers, flake, flake_label).animate.shift(UP*pick_up_height), run_time=1)

            held_layers.extend([flake, flake_label])
            stack_top_y += layer_spacing + line_thickness*0.01

        self.next_slide()
        self.play(VGroup(substrate, substrate_label).animate.shift(DOWN*2.5))
        bottom_y = substrate.get_top()[1] + layer_spacing
        stack_shift = held_layers[-1].get_bottom()[1] - bottom_y
        self.next_slide()
        self.play(Group(*held_layers).animate.shift(DOWN*stack_shift), run_time=1)

        self.next_slide()
        self.play(FadeOut(Group(polymer, polymer_label)), run_time=1)

        top_gate, top_gate_label = create_layer(gate_color, 0.7, "TG")
        top_gate.shift(UP*6)
        top_gate_label.shift(UP*6).set_color(BLACK)
        self.add(top_gate, top_gate_label)
        top_gate_x = top_gate_label.get_x()
        top_y_hover = held_layers[-2].get_top()[1] + layer_spacing
        self.next_slide()
        self.play(top_gate.animate.move_to([0,held_layers[0].get_y() + layer_spacing-0.1,0]),
                  top_gate_label.animate.move_to([top_gate_x,held_layers[0].get_top()[1] + layer_spacing,0]),
                  run_time=1)

        self.wait()

                # --- Gold contacts to TG, RG, and BG ---
        contact_color = GOLD
        contact_length = 5.3   # make it long so it looks like a wire
        contact_height = 0.2   # thin

        def make_contact(target_layer, label_text="Contact"):
            contact = Rectangle(
                width=contact_length,
                height=contact_height,
                fill_color=contact_color,
                fill_opacity=1,
                stroke_color=BLACK
            )
            y_pos = target_layer.get_y()
            # start offscreen left
            contact.move_to([-8, y_pos, 0])
            # label to the left of wire
            label = Text(label_text, font_size=20).next_to(contact, LEFT, buff=0.1).set_color(BLACK)
            return contact, label

        # Contacts for TG, RG, BG
        tg_contact, tg_label = make_contact(top_gate, "TG Contact")
        rg_contact, rg_label = make_contact(held_layers[4], "RG Contact")  # RG was 2nd flake
        bg_contact, bg_label = make_contact(held_layers[-2], "BG Contact") # BG was last flake before TG

        # Animate contacts sliding in from left and stopping at target
        self.next_slide()
        self.play(FadeIn(tg_contact, tg_label))
        self.play(tg_contact.animate.shift(RIGHT*3.5),)

        self.play(FadeIn(rg_contact, rg_label))
        self.play(rg_contact.animate.shift(RIGHT*2.5))

        self.play(FadeIn(bg_contact, bg_label))
        self.play(bg_contact.animate.shift(RIGHT*3.5))

        self.wait()

        visible = ImageMobject("visible.jpg").scale(0.45).move_to([-3,1.5,0])
        ir = ImageMobject("IR.jpg").scale(0.75).move_to([3,1.5,0])
        sMIM = ImageMobject("sMIM.jpg").scale(0.75).move_to([-3,1.5,0])
        cut = ImageMobject("Cut.jpg").scale(0.45).move_to([3,1.5,0])

        self.next_slide()
        self.play(FadeIn(visible))

        self.next_slide()
        self.play(FadeIn(ir))

        self.next_slide()
        self.play(FadeOut(visible), FadeIn(sMIM))

        self.next_slide()
        self.play(FadeOut(ir), FadeIn(cut))
        self.wait()

        held_layers.remove(polymer)
        held_layers.remove(polymer_label)
        stack = VGroup(held_layers, top_gate, top_gate_label, tg_contact, rg_contact, bg_contact, substrate, substrate_label)
        vector_magnet = SVGMobject("vector_magnet.svg").shift(LEFT*2.75+UP*1).scale(3)

        self.next_slide()
        self.play(stack.animate.scale(0.07).shift(RIGHT*0.6+UP*0.2), Write(vector_magnet), FadeOut(Group(sMIM, cut)))
        self.wait()

        mkelvin = MathTex(r"T\sim \text{mili-Kelvin}").shift(RIGHT*4).set_color(BLACK)
        self.next_slide()
        self.play(Write(mkelvin))
        self.wait()

        self.next_slide()
        self.play(*[FadeOut(mob) for mob in self.mobjects])
        self.wait()

        haoxin_transport = ImageMobject("Haoxin_PD.png")
        haoxin_cite = Tex("[5]").next_to(haoxin_transport,RIGHT).set_color(ORANGE).scale(0.8)
        self.play(FadeIn(haoxin_transport.scale(0.7)), FadeIn(haoxin_cite))
        self.wait()

        andrea = Group(
            ImageMobject("andrea_pic.png"),
            ImageMobject("andrea_name.png").scale(1.12)
        ).arrange(DOWN,buff=0).scale(1.05)

        people = Group(
            ImageMobject("ludwig.png"),
            ImageMobject("abhay.png"),
            ImageMobject("noah.png"),
            andrea
        ).arrange(RIGHT, buff=1).scale(0.3)#.shift(RIGHT*1.75 + DOWN*2)

        people_outlines = Group()
        for mob in people:
            people_outlines.add(SurroundingRectangle(mob, color=GRAY, buff=0))

        self.next_slide()
        self.play(*[FadeOut(mob) for mob in self.mobjects])
        self.play(FadeIn(people, people_outlines))
        self.wait()
        

        citations = VGroup(
            Text("[1]: 'Chiral superconductivity from repulsive interactions' in doped graphene — Rahul Nandkishore, Leonid Levitov, Andrey Chubukov (2011)").scale(0.3),
            Text("[2]: 'Anomalous Quantum Hall Effect: An Incompressible Quantum Fluid with Fractionally Charged Excitations'- R. B. Laughlin, Phys. Rev. Lett. 50, 1395 (1983)").scale(0.3),
            Text("[3]: 'Band structure of ABC-stacked graphene trilayers' - Fan Zhang, Bhagawan Sahu, Hongki Min, and A. H. MacDonald, Phys. Rev. B 82, 035409 (2010)").scale(0.3),
            Text("[4]: 'Electronic properties of multilayer graphene' -  M. Koshino and E. McCann, Phys. Rev. B 80, 165409 (2009)").scale(0.3),
            Text("[5]: 'Superconductivity in rhombohedral trilayer graphene' – Haoxin Zhou, Tian Xie, Takashi Taniguchi, Kenji Watanabe, Andrea F. Young (2021)").scale(0.3)
        ).arrange(DOWN, aligned_edge=LEFT).set_color(BLACK)

        self.next_slide()
        self.play(FadeOut(Group(people, people_outlines)))
        self.play(Write(citations))
        self.wait()



        



class Fabrication(Slide):
    def construct(self):

        self.camera.background_color = LIGHTER_GRAY
        # --- Colors ---
        polymer_color = GREY
        bn_color = BLUE
        graphene_color = PINK
        gate_color = PURPLE
        substrate_color = PURPLE_A

        # --- Parameters ---
        x_left, x_right = -3, 3
        line_thickness = 16
        layer_spacing = 0.25
        pick_up_height = 2

        # --- Substrate ---
        substrate = Line([-5,0,0],[5,0,0], stroke_width=line_thickness).set_color(substrate_color)
        substrate.shift(DOWN*0.5)
        substrate_label = Tex("SiO$_2$ Substrate", font_size=24).next_to(substrate, RIGHT, buff=.2).set_color(BLACK)
        self.add(substrate)
        self.add(substrate_label)

        # Track stack top
        stack_top_y = 0.1

        # --- Helper function for layers + labels ---
        def create_layer(color, width_scale=1.0, label_text=""):
            layer = Line([x_left*width_scale,0,0],[x_right*width_scale,0,0], stroke_width=line_thickness).set_color(color)
            label = Text(label_text, font_size=24).next_to(layer, RIGHT, buff=0.2)
            return layer, label

        # --- Step 0: Polymer above ---
        polymer, polymer_label = create_layer(polymer_color, 1.7, "Polymer")
        polymer.shift(UP*3)
        polymer_label.shift(UP*3).set_color(BLACK)
        self.add(polymer, polymer_label)
        held_layers = [polymer, polymer_label]

        self.wait(0.5)

        # --- Flake info: list of (color, label, width) ---
        flakes_info = [
            (bn_color, "hBN", 1.3),
            (graphene_color, "RG", 1),
            (bn_color, "hBN", 1.3),
            (gate_color, "BG", .7)
        ]

        # --- Pick up each flake ---
        for i, (color, label_text, width) in enumerate(flakes_info):

            self.next_slide()
            # Create flake
            flake, flake_label = create_layer(color, width, label_text)
            # Start offscreen right
            # flake.shift(RIGHT*8)
            flake_label.set_color(BLACK)
            # self.add(flake, flake_label)
            # # Fly in horizontally
            # self.play(flake.animate.shift(LEFT*8), flake_label.animate.shift(LEFT*8), run_time=1)

            self.play(FadeIn(VGroup(flake, flake_label)))

            # Polymer moves down just above current stack
            polymer_target_y = stack_top_y + layer_spacing
            self.play(Group(*held_layers).animate.shift(DOWN*(polymer.get_bottom()[1]-polymer_target_y)), run_time=1)

            # Lift polymer + held flakes together
            self.play(Group(*held_layers, flake, flake_label).animate.shift(UP*pick_up_height), run_time=1)

            # Add flake to held layers
            held_layers.extend([flake, flake_label])
            stack_top_y += layer_spacing + line_thickness*0.01

        # --- Drop entire stack onto substrate ---
        self.next_slide()
        self.play(VGroup(substrate, substrate_label).animate.shift(DOWN*2.5))
        # Compute shift to place bottom layer just above substrate
        bottom_y = substrate.get_top()[1] + layer_spacing
        stack_shift = held_layers[-1].get_bottom()[1] - bottom_y
        self.next_slide()
        self.play(Group(*held_layers).animate.shift(DOWN*stack_shift), run_time=1)

        # --- Fade out polymer only ---
        self.next_slide()
        self.play(FadeOut(Group(polymer, polymer_label)), run_time=1)

        # --- Top gate comes in from top ---
        top_gate, top_gate_label = create_layer(gate_color, 0.7, "TG")
        top_gate.shift(UP*6)
        top_gate_label.shift(UP*6).set_color(BLACK)
        self.add(top_gate, top_gate_label)
        top_gate_x = top_gate_label.get_x()
        # Move down to hover above top BN
        top_y_hover = held_layers[-2].get_top()[1] + layer_spacing
        self.next_slide()
        self.play(top_gate.animate.move_to([0,held_layers[0].get_y() + layer_spacing-0.1,0]),
                  top_gate_label.animate.move_to([top_gate_x,held_layers[0].get_top()[1] + layer_spacing,0]),
                  run_time=1)

        self.wait()

        visible = ImageMobject("visible.jpg").scale(0.45).move_to([-3,1.5,0])
        ir = ImageMobject("IR.jpg").scale(0.75).move_to([3,1.5,0])
        sMIM = ImageMobject("sMIM.jpg").scale(0.75).move_to([-3,1.5,0])
        cut = ImageMobject("Cut.jpg").scale(0.45).move_to([3,1.5,0])

        self.next_slide()
        self.play(FadeIn(visible))

        self.next_slide()
        self.play(FadeIn(ir))

        self.next_slide()
        self.play(FadeOut(visible), FadeIn(sMIM))

        self.next_slide()
        self.play(FadeOut(ir), FadeIn(cut))
        self.wait()
        
                # --- Gold contacts to TG, RG, and BG ---
        contact_color = GOLD
        contact_length = 5.3   # make it long so it looks like a wire
        contact_height = 0.2   # thin

        def make_contact(target_layer, label_text="Contact"):
            contact = Rectangle(
                width=contact_length,
                height=contact_height,
                fill_color=contact_color,
                fill_opacity=1,
                stroke_color=BLACK
            )
            y_pos = target_layer.get_y()
            # start offscreen left
            contact.move_to([-8, y_pos, 0])
            # label to the left of wire
            label = Text(label_text, font_size=20).next_to(contact, LEFT, buff=0.1).set_color(BLACK)
            return contact, label

        # Contacts for TG, RG, BG
        tg_contact, tg_label = make_contact(top_gate, "TG Contact")
        rg_contact, rg_label = make_contact(held_layers[4], "RG Contact")  # RG was 2nd flake
        bg_contact, bg_label = make_contact(held_layers[-2], "BG Contact") # BG was last flake before TG

        # Animate contacts sliding in from left and stopping at target
        self.next_slide()
        self.play(FadeIn(tg_contact, tg_label))
        self.play(tg_contact.animate.shift(RIGHT*3.5),)

        self.play(FadeIn(rg_contact, rg_label))
        self.play(rg_contact.animate.shift(RIGHT*2.5))

        self.play(FadeIn(bg_contact, bg_label))
        self.play(bg_contact.animate.shift(RIGHT*3.5))

        self.wait()

        held_layers.remove(polymer)
        held_layers.remove(polymer_label)

        stack = VGroup(held_layers, top_gate, top_gate_label, tg_contact, rg_contact, bg_contact, substrate, substrate_label)
        vector_magnet = SVGMobject("vector_magnet.svg").shift(LEFT*2.75+UP*1).scale(3)
        self.next_slide()
        self.play(stack.animate.scale(0.07).shift(RIGHT*0.6+UP*0.2), Write(vector_magnet), FadeOut(Group(sMIM, cut)))
        self.wait()


class CitationSlide(Slide):
    def construct(self):
        self.camera.background_color = LIGHTER_GRAY
        andrea = Group(
            ImageMobject("andrea_pic.png"),
            ImageMobject("andrea_name.png").scale(1.12)
        ).arrange(DOWN,buff=0).scale(1.05)

        people = Group(
            ImageMobject("ludwig.png"),
            ImageMobject("abhay.png"),
            ImageMobject("noah.png"),
            andrea
        ).arrange(RIGHT, buff=1).scale(0.3)#.shift(RIGHT*1.75 + DOWN*2)

        people_outlines = Group()
        for mob in people:
            people_outlines.add(SurroundingRectangle(mob, color=GRAY, buff=0))

        self.play(FadeIn(people, people_outlines))
        self.wait()
        

        citations = VGroup(
            Text("[1]: 'Chiral superconductivity from repulsive interactions' in doped graphene — Rahul Nandkishore, Leonid Levitov, Andrey Chubukov (2011)").scale(0.3),
            Text("[2]: 'Anomalous Quantum Hall Effect: An Incompressible Quantum Fluid with Fractionally Charged Excitations'- R. B. Laughlin, Phys. Rev. Lett. 50, 1395 (1983)").scale(0.3),
            Text("[3]: 'Band structure of ABC-stacked graphene trilayers' - Fan Zhang, Bhagawan Sahu, Hongki Min, and A. H. MacDonald, Phys. Rev. B 82, 035409 (2010)").scale(0.3),
            Text("[4]: 'Electronic properties of multilayer graphene' -  M. Koshino and E. McCann, Phys. Rev. B 80, 165409 (2009)").scale(0.3),
            Text("[5]: 'Superconductivity in rhombohedral trilayer graphene' – Haoxin Zhou, Tian Xie, Takashi Taniguchi, Kenji Watanabe, Andrea F. Young (2021)").scale(0.3)
        ).arrange(DOWN, aligned_edge=LEFT).set_color(BLACK)

        self.next_slide()
        self.play(FadeOut(Group(people, people_outlines)))
        self.play(Write(citations))
        self.wait()

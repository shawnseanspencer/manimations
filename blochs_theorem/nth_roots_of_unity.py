class NthRootofUnity(Scene):
    def construct(self):
        r = 2.5 #circle radius
        ghost_circle = Circle(radius=r).set_opacity(0)
        circular_crystal = VGroup()
        N = 12 #Number of points on circle
        for i in range(N):
            angle = 2*PI/N * i
            circular_crystal.add(
                Dot(color=BLUE).move_to(ghost_circle.get_center() + r*np.array([np.cos(angle), np.sin(angle), 0]))
            )
        return_circular_crystal = circular_crystal.copy() #What to finish this scene with
        self.add(circular_crystal)
        self.wait()
        plane = NumberPlane(
            x_range = [-3, 3],
            y_range =  [-3, 3],
            x_length = 6*r,
            y_length = 6*r
        ).set_opacity(0.3).set_color(GOLD)
        ax = Axes(
            x_range = [-2.8, 2.8],
            y_range =  [-1.6, 1.6],
            x_length = 5.6*r,
            y_length = 3.2*r
        ).set_opacity(0.4).set_color(WHITE)

        ax_labels = VGroup(
            Tex("Re").move_to([6.7,0.5,0]),
            Tex("Im").move_to([-0.8,3.7,0])
        )

        axes = VGroup(ax, plane, *[mob for mob in ax_labels])
        self.play(FadeIn(axes))
        self.wait(2)

        traveller = Dot(color=RED).scale(1.3)
        trail = TracedPath(traveller.get_center, dissipating_time=0.3, stroke_opacity=[0, 1], stroke_width=8).set_color(RED)
        self.add(traveller, trail)

        nth_roots = VGroup(
            MathTex("e^{i0\\pi/6"),
            MathTex("e^{i\\pi/6}"),
            MathTex("e^{i\\pi/3}"),
            MathTex("e^{i\\pi/2}"),
            MathTex("e^{i2\\pi/3}"),
            MathTex("e^{i5\\pi/6}"),
            MathTex("e^{i\\pi}"),
            MathTex("e^{i7\\pi/6}"),
            MathTex("e^{i4\\pi/3}"),
            MathTex("e^{i3\\pi/2}"),
            MathTex("e^{i5\\pi/3}"),
            MathTex("e^{i11\\pi/6}"),
        )


        angle = 2*PI/N
        for i, point in enumerate(circular_crystal):
            angle = 2*PI/N
            self.play(
                traveller.animate(path_arc=angle).move_to(point.get_center()),
                FadeIn(nth_roots[i].move_to(point.get_center()*1.2+RIGHT*0.2)),
                run_time=0.8
            )
        self.play(
            traveller.animate(path_arc=angle).move_to(circular_crystal[0].get_center())
        )
        self.wait()
        self.remove(trail)
        self.wait()
        nth_root_forms = VGroup(
            MathTex("e^{i2\\pi(0)/12}"),
            MathTex("e^{i2\\pi(1)/12}"),
            MathTex("e^{i2\\pi(2)/12}"),
            MathTex("e^{i2\\pi(3)/12}"),
            MathTex("e^{i2\\pi(4)/12}"),
            MathTex("e^{i2\\pi(5)/12}"),
            MathTex("e^{i2\\pi(6)/12}"),
            MathTex("e^{i2\\pi(7)/12}"),
            MathTex("e^{i2\\pi(8)/12}"),
            MathTex("e^{i2\\pi(9)/12}"),
            MathTex("e^{i2\\pi(10)/12}"),
            MathTex("e^{i2\\pi(11)/12}"),
        )
        for i, point in enumerate(circular_crystal):
            angle = 2*PI/N
            nth_root_forms[i].scale(0.8).move_to(point.get_center()*1.2+RIGHT*0.1)

        self.play(
            *[ReplacementTransform(nth_roots[i], nth_root_forms[i]) for i, root in enumerate(nth_roots)]
        )
        self.wait(2)
        nth_root_statement = Tex("""Nth Roots\n
                                 Of Unity""")
        nth_root_box = BackgroundRectangle(nth_root_statement)
        nth_root_rect = SurroundingRectangle(nth_root_statement).set_color(GOLD)
        nth_root_statement_group = VGroup(nth_root_box, nth_root_statement, nth_root_rect)
        self.play(FadeIn(nth_root_statement_group))
        self.wait()

        nth_root_form = MathTex("e^{i2\\pi s/N}").shift(DOWN*0.55)
        nth_root_form_box = BackgroundRectangle(nth_root_form).set_opacity(0.8)
        nth_root_form_group = VGroup(nth_root_form_box, nth_root_form)
        self.play(
            FadeIn(nth_root_form_group),
            nth_root_statement_group.animate.shift(UP*0.3),
            ReplacementTransform(nth_root_rect, SurroundingRectangle(VGroup(nth_root_statement.copy().shift(UP*.3), nth_root_form)).set_color(GOLD)),
            ReplacementTransform(nth_root_box, BackgroundRectangle(VGroup(nth_root_statement.copy().shift(UP*.3), nth_root_form)).set_color(BLACK))
        )
        self.wait()

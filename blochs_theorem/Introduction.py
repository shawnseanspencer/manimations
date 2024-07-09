class Introduction(ThreeDScene):
    def construct(self):
        lattice_points_1D = VGroup()
        ah_lattice = 1 #Lattice constant horizontal
        av_lattice = 1 #Lattice constant vertical
        point_count_h = 16 #Horizontal point count

        for i in range(point_count_h):
            lattice_points_1D.add(Dot(color=BLUE).move_to(LEFT*10+RIGHT*ah_lattice*i))
        lattice_points_1D = lattice_points_1D.move_to([0,0,0])

        self.play(Write(lattice_points_1D), run_time=2)
        infinite_arrows = VGroup(
            Arrow(start=[-1,0,0], end=[-4,0,0]),
            Arrow(start=[1,0,0], end=[4,0,0])
        ).shift(UP*1)
        infinite_arrows.add(MathTex("\\infty").next_to(infinite_arrows[0].get_end(), LEFT))
        infinite_arrows.add(MathTex("\\infty").next_to(infinite_arrows[1].get_end(), RIGHT))

        self.wait() 
        self.play(Write(infinite_arrows))
        self.wait()
        self.play(lattice_points_1D.animate.shift(RIGHT*1.5*ah_lattice), run_time=3)
        self.wait()
        self.play(lattice_points_1D.animate.shift(RIGHT*1.5*ah_lattice), run_time=1.5)
        self.play(Unwrite(infinite_arrows))
        self.wait()

        point_count_v = 17 #Vertical point coint
        lattice_points_2D = VGroup()
        for i in range(point_count_v):
            lattice_points_2D.add(
                lattice_points_1D.copy().shift(UP*4.5+DOWN*av_lattice*i)
            )
        
        lattice_points_2D.move_to([0,0,0])
        self.play(Write(lattice_points_2D), FadeOut(lattice_points_1D), run_time=4)
        self.remove(lattice_points_1D)
        self.wait()

        lattice_constant_val = ValueTracker(np.linalg.norm(lattice_points_2D[0][0].get_center() - lattice_points_2D[0][1].get_center())) 
        lattice_constant = VGroup(
            MathTex("a"),
            DecimalNumber(lattice_constant_val.get_value(), num_decimal_places=3, include_sign=False, unit=None)
        ).move_to([4.2,2,0]).set_color(RED).scale(1.5)
        lattice_constant[1] = lattice_constant[1].next_to(lattice_constant[0], DOWN, buff=0.2)

        lattice_constant[1].add_updater(
            lambda d: d.set_value(lattice_constant_val.get_value())
        )

        lattice_bracket = BraceBetweenPoints(lattice_points_2D[10][10].get_center(), lattice_points_2D[10][11].get_center()).set_color(RED)
    
        self.play(Write(lattice_constant), DrawBorderThenFill(lattice_bracket))

        def track_lattice_constant(scale_factor):
            return (lattice_constant_val.animate.set_value(np.linalg.norm(lattice_points_2D.copy().scale(scale_factor)[0][0].get_center() 
                                                                          - lattice_points_2D.copy().scale(scale_factor)[0][1].get_center())))

        self.wait(2)
        self.play(FadeOut(lattice_bracket))
        self.play(
            lattice_points_2D.animate.scale(0.9), 
            track_lattice_constant(0.9),
            run_time=1.5,
        )
        self.wait()
        self.play(Rotate(lattice_points_2D, about_point=ORIGIN), run_time=5)
        self.play(Unwrite(lattice_constant))
        self.wait()
        self.move_camera(
            phi = -100*DEGREES,
            theta = 70*DEGREES,
            zoom = 0.25,
            frame_center = [0,0,-5],
            run_time=2
        )
        
        self.wait()
        threeD_lattice_layer = VGroup()
        for i in range(len(lattice_points_2D)):
            for j in range(len(lattice_points_2D[0])):
                threeD_lattice_layer.add(Dot().move_to(lattice_points_2D[i][j].get_center()).scale(1).set_color(BLUE))
                
        
        self.play(Create(threeD_lattice_layer))

        threeD_lattice = VGroup()
        for i in range(12):
            threeD_lattice.add(threeD_lattice_layer.copy().shift(np.array([0,0,-1])*i*lattice_constant_val.get_value()))

        self.wait()

        self.play(Create(threeD_lattice))

        self.wait()
        self.move_camera(
            phi = -120*DEGREES,
            run_time = 6
        )

        self.wait()

        self.move_camera(
            zoom = 4
        )

        self.wait()
        self.move_camera(
            theta = 50*DEGREES,
            phi = -100*DEGREES,
            rum_time = 7
        )
        self.wait()

        ghost_lattice = Group(*[mob.copy().set_opacity(0.2) for mob in self.mobjects])
        not_ghost = Group(*[mob for mob in self.mobjects])

        self.add(ghost_lattice)
        self.remove(not_ghost)
        self.play(
            *[mob.animate.shift(([0, 0, lattice_constant_val.get_value()])) for mob in not_ghost],
            run_time = 5
        )
        self.play(
            *[mob.animate.shift(([0 ,-lattice_constant_val.get_value(), 0])) for mob in not_ghost],
            run_time = 5
        )
        self.wait()
        #Pointing out finiteness of lattice
        cam_phi, cam_theta = self.camera.phi, self.camera.theta
        cam_frame_center = self.camera.frame_center
        relative_position = cam_frame_center -  np.array([
            np.sin(cam_phi) * np.cos(cam_theta),
            np.sin(cam_phi) * np.sin(cam_theta),
            np.cos(cam_phi)
        ]) 

        self.wait()
        self.remove(ghost_lattice)

        self.wait()

        self.move_camera(
            frame_center = [0,0,0],
            phi = 0*DEGREES,
            theta = -90*DEGREES,
            zoom = 1,
            run_time = 5
        )
        self.wait()
        self.play(
            FadeOut(threeD_lattice),
        )
        self.play(FadeIn(lattice_points_1D.move_to([0,0,0])), FadeOut(lattice_points_2D), FadeOut(threeD_lattice_layer))
        self.wait()

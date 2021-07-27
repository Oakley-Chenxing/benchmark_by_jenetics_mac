package io.jenetics.griewank;

import io.jenetics.Mutator;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.engine.Limits;
import io.jenetics.util.ISeq;
import io.jenetics.util.RandomRegistry;

import io.jenetics.ext.SingleNodeCrossover;
import io.jenetics.ext.util.TreeNode;

import io.jenetics.prog.ProgramGene;
import io.jenetics.prog.op.EphemeralConst;
import io.jenetics.prog.op.MathExpr;
import io.jenetics.prog.op.MathOp;
import io.jenetics.prog.op.Op;
import io.jenetics.prog.op.Var;
import io.jenetics.prog.regression.Error;
import io.jenetics.prog.regression.LossFunction;
import io.jenetics.prog.regression.Regression;
import io.jenetics.prog.regression.Sample;
import java.util.List;
import java.util.function.Function;
import static java.util.stream.Collectors.toList;

import io.jenetics.Chromosome;
import io.jenetics.Genotype;
import io.jenetics.Mutator;
import io.jenetics.engine.Codec;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.util.ISeq;
import io.jenetics.util.RandomRegistry;

import io.jenetics.ext.SingleNodeCrossover;

import io.jenetics.prog.ProgramChromosome;

public class Griewank6 {
	public static final double XMIN = -10.0;
	public static final double XMAX =  10.0;
	public static final Op<Double> PDIV = Op.of("pdiv", 2, v -> (v[1] == 0) ? 1 : v[0]/v[1]);

	// Definition of the allowed operations.
	static final ISeq<Op<Double>> OPERATIONS = ISeq.of(
		MathOp.ADD,
		MathOp.SUB,
		MathOp.MUL,
		PDIV
	);

	// Definition of the terminals.
	static final ISeq<Op<Double>> TERMINALS = ISeq.of(
		EphemeralConst.of(() -> RandomRegistry.random().nextDouble())
	);

	static final Codec<ISeq<Function<Double[], Double>>, ProgramGene<Double>> CODEC =
		Codec.of(
			Genotype.of(
				// First 'program'
				ProgramChromosome.of(
					5,
					ch -> ch.root().size() <= 50,
					OPERATIONS,
					TERMINALS
				),
				// Second 'program'
				ProgramChromosome.of(
					5,
					ch -> ch.root().size() <= 50,
					OPERATIONS,
					TERMINALS
				),
	            // Third 'program'
	            ProgramChromosome.of(
	                    5,
	                    ch -> ch.root().size() <= 50,
	                    OPERATIONS,
	                    TERMINALS
	            ),
	            // 4th 'program'
	            ProgramChromosome.of(
	                   5,
	                    ch -> ch.root().size() <= 50,
	                    OPERATIONS,
	                    TERMINALS
	            ),
	            // 5th 'program'
	            ProgramChromosome.of(
	                    5,
	                    ch -> ch.root().size() <= 50,
	                    OPERATIONS,
	                    TERMINALS
	            ),
	            // 6th 'program'
	            ProgramChromosome.of(
	                    5,
	                    ch -> ch.root().size() <= 50,
	                    OPERATIONS,
	                    TERMINALS
	            )
			),
			gt -> gt.stream()
				.map(c -> (ProgramGene<Double>) c.gene())
				.collect(ISeq.toISeq())
		);

	static double fitness(final ISeq<Function<Double[], Double>> programs) {
		assert programs.size() == 6;
		double[] rep = programs.stream()
				.map(p -> p.apply(new Double[0]))
				.map(x -> Math.max(x, XMIN))
				.map(x -> Math.min(x, XMAX))
				.mapToDouble(d -> d).toArray();
		
		double fit = 0.0;
		double sum = 0.0;
		double prod = 1.0;
		for (int t=0; t<rep.length; t++) {
			double temp1 = rep[t] * rep[t];
			sum += temp1;
			double temp2 = (double)(Math.cos(rep[t]/Math.sqrt(t+1)));
			prod *= temp2;
			fit += (1/4000)*sum - prod + 1.0;
		}

		return fit;
	}
	
//	 float sum  = 0.0f;
//	 float prod = 1.0f;
//	   for (int i=0; i<rep.length; i++) {
//	       System.out.println("rep["+i+"]: "+rep[i]);
//		   float temp1 = rep[i]*rep[i];
//		   System.out.println("temp1: "+temp1);
//		   sum += temp1;
//		   System.out.println("sum: "+sum);
//		   float temp2 = (float) (Math.cos(rep[i]/Math.sqrt(i+1)));
//		   System.out.println("temp2: "+temp2);
//		   prod *= temp2;
//		   System.out.println("prod: "+prod);
//}
//		return (1/4000)*sum - prod + 1.0f;
//	
	

	public static void main(final String[] args) {
		long startTime = System.currentTimeMillis();
		final Engine<ProgramGene<Double>, Double> engine = Engine
			.builder(Griewank6::fitness, CODEC)
			.minimizing()
			.alterers(
				new SingleNodeCrossover<>(),
				new Mutator<>())
			.build();

		final Genotype<ProgramGene<Double>> gt = engine.stream()
			.limit(20000)
			.collect(EvolutionResult.toBestGenotype());

		final ISeq<Function<Double[], Double>> programs = CODEC.decode(gt);
		
		List<Double> best = programs.stream()
									.map(p -> p.apply(new Double[0]))
									.map(x -> Math.max(x, XMIN))
									.map(x -> Math.min(x, XMAX))
									.collect(toList());
		
		System.out.format("%.6f\n", fitness(programs));
		System.out.println(best.stream().map(x -> String.format("%.6f", x)).collect(toList()));
		long endTime = System.currentTimeMillis();
        System.out.println(endTime-startTime + " ms");
	}
}

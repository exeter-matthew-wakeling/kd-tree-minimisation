import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Reads a file containing a set of points each containing some number of
 * "in" values and "out" values, and predicts the "out" values from the
 * "in" values. The points will be placed in a k-d tree. Predictions can
 * be made by finding the n most similar points and taking the average of
 * the "out" values.
 *
 * @author Matthew Wakeling
 */
public class PredictFromStats
{
	/**
	 * Main program.
	 *
	 * @param args the command-line arguments.
	 */
	public static void main(String args[]) throws IOException {
		List<Integer> columns = new ArrayList<Integer>();
		List<Limits> range = new ArrayList<Limits>();
		List<Integer> valueColumns = new ArrayList<Integer>();
		int pop = 10;
		int dbLimit = Integer.MAX_VALUE;
		double maxradius = 5.0;
		String lookupFile = null;
		boolean printRetrieved = false;
		for (int i = 0; i < args.length; i++) {
			if (args[i].startsWith("-v")) {
				int columnCount = Integer.parseInt(args[i].substring(2));
				for (int o = 0; o < columnCount; o++) {
					int column = Integer.parseInt(args[o + i + 1]);
					valueColumns.add(column);
					range.add(new Limits());
				}
				i += columnCount;
			} else if (args[i].startsWith("-c")) {
				int columnCount = Integer.parseInt(args[i].substring(2));
				for (int o = 0; o < columnCount; o++) {
					int column = Integer.parseInt(args[o + i + 1]);
					columns.add(column);
					range.add(new Limits());
				}
				i += columnCount;
			} else if (args[i].startsWith("-l")) {
				i += 2;
			} else if ("--limit".equals(args[i])) {
				dbLimit = Integer.parseInt(args[i + 1]);
				i++;
			} else if ("-p".equals(args[i])) {
				pop = Integer.parseInt(args[i + 1]);
				i++;
			} else if ("-b".equals(args[i])) {
				maxradius = Double.parseDouble(args[i + 1]) / 1000.0;
				i++;
			} else if ("--lookup".equals(args[i])) {
				lookupFile = args[i + 1];
				i++;
			} else if ("--retrieved".equals(args[i])) {
				printRetrieved = true;
			} else {
				throw new IllegalArgumentException("Invalid argument " + args[i]);
			}
		}

		List<Point> points = new ArrayList<Point>();
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		String row = in.readLine();
		while ((row != null) && (points.size() < dbLimit)) {
			List<Double> vs = new ArrayList<Double>();
			for (int i = 0; i < valueColumns.size(); i++) {
				vs.add(Double.NaN);
			}
			List<Double> ins = new ArrayList<Double>();
			for (int i = 0; i < columns.size(); i++) {
				ins.add(Double.NaN);
			}
			row = row.trim() + " ";
			int val = 1;
			while (!"".equals(row.trim())) {
				int spacePos = row.indexOf(' ');
				int tabPos = row.indexOf('\t');
				if ((tabPos > -1) && (tabPos < spacePos)) {
					spacePos = tabPos;
				}
				for (int i = 0; i < valueColumns.size(); i++) {
					if (val == valueColumns.get(i)) vs.set(i, new Double(row.substring(0, spacePos)));
				}
				for (int i = 0; i < columns.size(); i++) {
					if (val == columns.get(i)) ins.set(i, new Double(row.substring(0, spacePos)));
				}
				row = row.substring(spacePos).trim() + " ";
				val++;
			}
			boolean good = true;
			for (int i = 0; i < valueColumns.size(); i++) {
				good = good && (!Double.isNaN(vs.get(i)));
			}
			for (int i = 0; i < columns.size(); i++) {
				good = good && (!Double.isNaN(ins.get(i)));
			}
			if (good) {
				points.add(new Point(vs, ins));
				for (int i = 0; i < columns.size(); i++) {
					range.get(i).expand(ins.get(i));
				}
			}
			row = in.readLine();
		}

		for (int i = 0; i < args.length; i++) {
			if ("-v".equals(args[i])) {
				i++;
			} else if (args[i].startsWith("-c")) {
				int columnCount = Integer.parseInt(args[i].substring(2));
				i += columnCount;
			} else if (args[i].startsWith("-l")) {
				int columnNo = Integer.parseInt(args[i].substring(2));
				range.set(columnNo + 1, new Limits(Double.parseDouble(args[i + 1]), Double.parseDouble(args[i + 2])));
				i += 2;
			} else if ("-p".equals(args[i])) {
				i++;
			} else if ("-b".equals(args[i])) {
				i++;
			}
		}

/*		for (Point p : points) {
			List<Double> ins = p.getIns();
			for (int i = 0; i < ins.size(); i++) {
				ins.set(i, (ins.get(i) - range.get(i).getMin()) / (range.get(i).getMax() - range.get(i).getMin()));
			}
		}*/

		Tree tree = new Tree(new ArrayList<Point>(points));

		if (lookupFile != null) {
			points = new ArrayList<Point>();
			in = new BufferedReader(new FileReader(lookupFile));
			row = in.readLine();
			while (row != null) {
				List<Double> vs = new ArrayList<Double>();
				for (int i = 0; i < valueColumns.size(); i++) {
					vs.add(Double.NaN);
				}
				List<Double> ins = new ArrayList<Double>();
				for (int i = 0; i < columns.size(); i++) {
					ins.add(Double.NaN);
				}
				row = row.trim() + " ";
				int val = 1;
				while (!"".equals(row.trim())) {
					int spacePos = row.indexOf(' ');
					int tabPos = row.indexOf('\t');
					if ((tabPos > -1) && (tabPos < spacePos)) {
						spacePos = tabPos;
					}
					for (int i = 0; i < valueColumns.size(); i++) {
						if (val == valueColumns.get(i)) vs.set(i, new Double(row.substring(0, spacePos)));
					}
					for (int i = 0; i < columns.size(); i++) {
						if (val == columns.get(i)) {
							//ins.set(i, (Double.parseDouble(row.substring(0, spacePos)) - range.get(i).getMin()) / (range.get(i).getMax() - range.get(i).getMin()));
							ins.set(i, Double.parseDouble(row.substring(0, spacePos)));
						}
					}
					row = row.substring(spacePos).trim() + " ";
					val++;
				}
				boolean good = true;
				for (int i = 0; i < valueColumns.size(); i++) {
					good = good && (!Double.isNaN(vs.get(i)));
				}
				for (int i = 0; i < columns.size(); i++) {
					good = good && (!Double.isNaN(ins.get(i)));
				}
				if (good) {
					points.add(new Point(vs, ins));
				}
				row = in.readLine();
			}
		}
		ParallelLookup lookup = new ParallelLookup(tree, pop, 0.0, maxradius, points, printRetrieved);
		while (lookup.hasNext()) {
			Lookup l = lookup.next();
			Point p = l.getIn();
			List<Stats> stats = l.getOut();
			if (!printRetrieved) {
				if (stats == null) {
					for (int i = 0; i < p.getVs().size(); i++) {
						if (i != 0) {
							System.out.print("\t");
						}
						System.out.print(p.getVs().get(i) + "\t0.0\t0.0");
					}
				} else {
					for (int i = 0; i < p.getVs().size(); i++) {
						if (i != 0) {
							System.out.print("\t");
						}
						System.out.print(p.getVs().get(i) + "\t" + stats.get(i).getMean() + "\t" + stats.get(i).getStdDev());
					}
				}
				System.out.println("");
			}
		}
	}

	public static class Lookup
	{
		private Point in;
		private List<Stats> out;

		public Lookup(Point in, List<Stats> out) {
			this.in = in;
			this.out = out;
		}

		public Point getIn() {
			return in;
		}

		public List<Stats> getOut() {
			return out;
		}
	}

	public static class ParallelLookup implements Iterator<Lookup>
	{
		private Tree tree;
		private int pop, outPos, inPos, end, aheadLimit;
		private double minradius, maxradius;
		private List<Point> inList;
		private TreeMap<Integer, Lookup> buffer = new TreeMap<Integer, Lookup>();
		boolean printRetrieved;

		public ParallelLookup(Tree tree, int pop, double minradius, double maxradius, List<Point> inList, boolean printRetrieved) {
			this.tree = tree;
			this.pop = pop;
			this.outPos = 0;
			this.inPos = 0;
			this.minradius = minradius;
			this.maxradius = maxradius;
			this.inList = inList;
			this.end = inList.size();
			this.aheadLimit = 4 * Runtime.getRuntime().availableProcessors();
			this.printRetrieved = printRetrieved;

			for (int i = 0; i < (printRetrieved ? 1 : Runtime.getRuntime().availableProcessors()); i++) {
				Thread t = new Thread(new ParallelLookupThread());
				t.start();
			}
		}

		public boolean hasNext() {
			return outPos < end;
		}

		public synchronized Lookup next() {
			Lookup retval = buffer.get(outPos);
			while (retval == null) {
				try {
					wait();
				} catch (InterruptedException e) {
				}
				retval = buffer.get(outPos);
			}
			outPos++;
			notifyAll();
			return retval;
		}

		public void remove() {
			throw new UnsupportedOperationException();
		}

		public synchronized int getNextNo() {
			while (inPos > outPos + aheadLimit) {
				try {
					wait();
				} catch (InterruptedException e) {
				}
			}
			return inPos++;
		}

		public synchronized void putResult(int nextNo, Lookup lookup) {
			buffer.put(nextNo, lookup);
			notifyAll();
		}

		public class ParallelLookupThread implements Runnable
		{
			public void run() {
				int nextNo = getNextNo();
				while (nextNo < end) {
					Point p = inList.get(nextNo);
					List<Stats> stats = tree.getStatsFor(p.getIns(), pop, minradius, maxradius, printRetrieved);
					putResult(nextNo, new Lookup(p, stats));
					nextNo = getNextNo();
				}
			}
		}
	}

	public static abstract class Node
	{
		public abstract Point getTrivialNearest(List<Double> centre);
		public abstract List<Point> getInRectangle(List<Limits> limits);
		public abstract double distFrom(List<Double> d);
	}

	public static class Tree extends Node
	{
		private List<Limits> limits;
		private double pivot;
		private int direction;
		private Node low, high;

		public Tree(List<Point> buckets) {
			this(0, buckets);
		}

		public Tree(int direction, List<Point> buckets) {
			this.direction = direction;
			build(buckets);
			limits = new ArrayList<Limits>();
			int dimensions = buckets.get(0).getIns().size();
			for (int i = 0; i < dimensions; i++) {
				limits.add(new Limits());
			}
			for (Point p : buckets) {
				for (int i = 0; i < dimensions; i++) {
					limits.get(i).expand(p.getIns().get(i));
				}
			}
		}

		public void build(List<Point> buckets) {
			int dimensions = buckets.get(0).getIns().size();
			List<Point> subLow = new ArrayList<Point>();
			List<Point> subHigh = new ArrayList<Point>();
			Collections.sort(buckets, new Comparator<Point>() {
						public int compare(Point a, Point b) {
							if (a.getIns().get(direction) > b.getIns().get(direction)) {
							return 1;
							} else if (a.getIns().get(direction) < b.getIns().get(direction)) {
							return -1;
							} else {
							return 0;
							}
						}
					}
					);
			int pos = buckets.size() / 2;
			pivot = buckets.get(pos).getIns().get(direction);
			/*if (pivot == buckets.get(0).getX()) pivot = Math.nextAfter(pivot, pivot + 1.0);
			for (Point b : buckets) {
				(b.getX() < pivot ? subLow : subHigh).add(b);
			}*/
			subLow.addAll(buckets.subList(0, pos));
			subHigh.addAll(buckets.subList(pos, buckets.size()));
			if (subLow.isEmpty() || subHigh.isEmpty()) {
				System.err.println("Zero-sized split, direction " + direction + ", low " + subLow.size() + ", high " + subHigh.size() + ", low " + subLow + ", high " + subHigh);
				direction = (direction + 1) % dimensions;
				build(buckets);
			} else {
				if (subLow.size() == 1) {
					low = subLow.get(0);
				} else {
					low = new Tree((direction + 1) % dimensions, subLow);
				}
				if (subHigh.size() == 1) {
					high = subHigh.get(0);
				} else {
					high = new Tree((direction + 1) % dimensions, subHigh);
				}
			}
		}

		public Point getTrivialNearest(List<Double> d) {
			return (d.get(direction) < pivot ? low : high).getTrivialNearest(d);
		}

		public Point getNearest(List<Double> d) {
			Point trivial = getTrivialNearest(d);
			double dist = trivial.distFrom(d);
			List<Limits> rectangle = new ArrayList<Limits>();
			for (int i = 0; i < d.size(); i++) {
				rectangle.add(new Limits(d.get(i) - dist, d.get(i) + dist));
			}
			List<Point> few = getInRectangle(rectangle);
			Point nearest = null;
			double nDist = Double.MAX_VALUE;
			for (Point b : few) {
				dist = b.distFrom(d);
				if (dist < nDist) {
					nDist = dist;
					nearest = b;
				}
			}
			return nearest;
		}

		public List<Point> getInRectangle(List<Limits> rectangle) {
			if (rectangle.get(direction).getMax() < pivot) {
				return low.getInRectangle(rectangle);
			} else if (rectangle.get(direction).getMin() >= pivot) {
				return high.getInRectangle(rectangle);
			} else {
				List<Point> retval = new ArrayList<Point>();
				retval.addAll(low.getInRectangle(rectangle));
				retval.addAll(high.getInRectangle(rectangle));
				return retval;
			}
		}

		public double distFrom(List<Double> d) {
			double sum = 0.0;
			for (int i = 0; i < d.size(); i++) {
				double da = d.get(i);
				Limits l = limits.get(i);
				if (da < l.getMin()) {
					sum += (da - l.getMin()) * (da - l.getMin());
				} else if (da > l.getMax()) {
					sum += (da - l.getMax()) * (da - l.getMax());
				}
			}
			return Math.sqrt(sum);
		}

		/**
		 * Returns the statistics for a given position. Note that the coordinates are in normalised
		 * space. Implements a nearest neighbour queue, finding usually pop nearest neighbours,
		 * but at least those within minradius and definitely not those outside maxradius.
		 *
		 * @param d the coordinates
		 * @param pop all else considered, perform stats on this many nearest neighbours
		 * @param minradius include all points within this radius
		 * @param maxradius disclude all points outside this radius
		 * @return a Stats object
		 */
		public List<Stats> getStatsFor(final List<Double> d, int pop, double minradius, double maxradius, boolean printRetrieved) {
			TreeSet<Node> queue = new TreeSet<Node>(new Comparator<Node>() {
					public int compare(Node a, Node b) {
					double ad = a.distFrom(d);
					double bd = b.distFrom(d);
					if (ad > bd) {
					return 1;
					} else if (ad < bd) {
					return -1;
					} else if (a == b) {
					return 0;
					} else {
					return a.hashCode() - b.hashCode();
					}
					}
					});
			queue.add(this);
			boolean notFinished = true;
			List<Point> points = new ArrayList<Point>();
			while (notFinished) {
				Node head = queue.pollFirst();
				if (head instanceof Point) {
					double dist = head.distFrom(d);
					if ((dist <= maxradius) && ((points.size() < pop) || (dist <= minradius))) {
						points.add((Point) head);
					} else {
						notFinished = false;
					}
				} else if (head instanceof Tree) {
					queue.add(((Tree) head).low);
					queue.add(((Tree) head).high);
				} else {
					notFinished = false;
				}
			}
			if (!points.isEmpty()) {
				int origPointCount = points.size();
				/*{
					double ssum = 0.0;
					int closestP = -1;
					double closest = Double.MAX_VALUE;
					int furthestP = -1;
					double furthest = -1.0;
					double[] sum = null;
					for (Point p : points) {
						List<Double> vs = p.getVs();
						if (sum == null) {
							sum = new double[vs.size()];
						}
						for (int i = 0; i < vs.size(); i++) {
							sum[i] += vs.get(i);
						}
					}
					for (int i = 0; i < sum.length; i++) {
						sum[i] = sum[i] / points.size();
					}
					for (int pointNo = 0; pointNo < points.size(); pointNo++) {
						double pSsum = 0.0;
						Point p = points.get(pointNo);
						List<Double> vs = p.getVs();
						for (int i = 0; i < vs.size(); i++) {
							double diff = vs.get(i) - sum[i];
							pSsum += diff * diff;
						}
						if (pSsum > furthest) {
							furthest = pSsum;
							furthestP = pointNo;
						}
						if (pSsum < closest) {
							closest = pSsum;
							closestP = pointNo;
						}
						ssum += pSsum;
					}
					ssum = ssum / points.size();
					//System.err.println("Closest = " + closest + ", ssum = " + ssum + ", sum.length = " + sum.length + ", points.size() = " + points.size() + ", furthest = " + furthest);
					if (closest > ssum * sum.length / points.size()) {
						// Use just the closest point.
						if (printRetrieved) {
							for (int pointNo = 1; pointNo < points.size(); pointNo++) {
								Point p = points.get(pointNo);
								System.out.print("0");
								List<Double> vs = p.getVs();
								for (int i = 0; i < vs.size(); i++) {
									System.out.print("\t" + vs.get(i));
								}
								for (int i = 0; i < sum.length; i++) {
									System.out.print("\t" + sum[i]);
								}
								System.out.println("\t" + (points.size() + 1));
							}
						}
						Point p = points.get(0);
						points = new ArrayList<Point>();
						points.add(p);
					}
				}*/
				//System.err.println("Points reduced from " + origPointCount + " to " + points.size());
				double[] sum = null;
				double[] ssum = null;
				double[] count = null;
				for (Point p : points) {
					List<Double> vs = p.getVs();
					if (printRetrieved) {
						System.out.print("1");
						for (int i = 0; i < vs.size(); i++) {
							System.out.print("\t" + vs.get(i));
						}
						System.out.println("");
					}
					if (sum == null) {
						sum = new double[vs.size()];
						ssum = new double[vs.size()];
						count = new double[vs.size()];
					}
					for (int i = 0; i < vs.size(); i++) {
						sum[i] += vs.get(i);
						ssum[i] += vs.get(i) * vs.get(i);
						count[i] += 1.0;
					}
				}
				Point furthest = points.get(points.size() - 1);
				double dist = furthest.distFrom(d);
				List<Stats> retval = new ArrayList<Stats>();
				for (int i = 0; i < sum.length; i++) {
					if (ssum[i] < (sum[i] * sum[i] / count[i])) {
						//System.err.println("ssum = " + ssum + ", sum * sum / count = " + (sum * sum / count));
						retval.add(new Stats(sum[i] / count[i], 0.0, Math.PI * dist * dist, count[i]));
					} else {
						retval.add(new Stats(sum[i] / count[i], Math.sqrt((ssum[i] - (sum[i] * sum[i]) / count[i]) / count[i]), Math.PI * dist * dist, count[i]));
					}
				}
				return retval;
			} else {
				return null;
			}
		}

		public String toString() {
			return "Direction " + direction + ", pivot: " + pivot + " (" + limits + ")";
		}
	}

	public static class Point extends Node
	{
		private List<Double> vs;
		private List<Double> ins;

		public Point(List<Double> vs, List<Double> ins) {
			this.vs = vs;
			this.ins = ins;
		}

		public List<Double> getIns() {
			return ins;
		}

		public List<Double> getVs() {
			return vs;
		}

		public Point getTrivialNearest(List<Double> centre) {
			return this;
		}

		public List<Point> getInRectangle(List<Limits> limits) {
			boolean in = true;
			for (int i = 0; i < limits.size(); i++) {
				in = in && limits.get(i).contains(ins.get(i));
			}
			if (in) {
				return Collections.singletonList(this);
			} else {
				return Collections.emptyList();
			}
		}

		public double distFrom(List<Double> d) {
			double sum = 0.0;
			for (int i = 0; i < d.size(); i++) {
				sum += (ins.get(i) - d.get(i)) * (ins.get(i) - d.get(i));
			}
			return Math.sqrt(sum);
		}

		public String toString() {
			StringBuilder retval = new StringBuilder("(");
			boolean needComma = false;
			for (double d : ins) {
				if (needComma) {
					retval.append(", ");
				}
				needComma = true;
				retval.append(d);
			}
			retval.append(")");
			return retval.toString();
		}
	}

	public static class Limits
	{
		private double min, max;

		public Limits() {
			this.min = Double.MAX_VALUE;
			this.max = -Double.MAX_VALUE;
		}

		public Limits(double min, double max) {
			this.min = min;
			this.max = max;
		}

		public void expand(double v) {
			min = Math.min(min, v);
			max = Math.max(max, v);
		}

		public double getMin() {
			return min;
		}

		public double getMax() {
			return max;
		}

		public boolean contains(double d) {
			return (d >= min) && (d < max);
		}

		public String toString() {
			return min + " - " + max;
		}
	}

	public static class Stats
	{
		private double mean, stdDev, area, count;

		public Stats(double mean, double stdDev, double area, double count) {
			this.mean = mean;
			this.stdDev = stdDev;
			this.area = area;
			this.count = count;
		}

		public double getMean() {
			return mean;
		}

		public double getStdDev() {
			return stdDev;
		}

		public double getArea() {
			return area;
		}

		public double getCount() {
			return count;
		}

		public String toString() {
			return "Mean: " + mean + ", Stddev: " + stdDev + ", area: " + area + ", count: " + count;
		}
	}
}

package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Lists;

import java.util.AbstractCollection;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import static com.google.common.base.Preconditions.checkElementIndex;

/**
 * Implementation of a persistent list. I.e 'set(i, val)' returns a new PersistentArrayList,
 * leaving the original PersistentArrayList un(perceptably)changed.
 *
 * Based on "A Persistent Union-Find Data Structure" Conchon & Filliaˆtre, 2007
 *
 * NOT THREADSAFE.
 * Accesses to a PersistentArrayList can potentially mutate a chain of PersistentArrayLists
 * in a thread-unsafe way during the method call.
 * TODO: add locks?
 *
 * Not resizable.
 *
 * @author sthomson@cs.cmu.edu
 */
public class PersistentArrayList<T> extends AbstractCollection<T> {
	// either an array or a diff chain based on an array
	private Data<T> data;
	private int size;

	/** Private class to hold our data. Either an array, or a diff chain based on an array */
	private static interface Data<T> {
		public T get(int i);
	}
	/** Just a plain old ArrayList */
	private static class ArrayData<T> implements Data<T> {
		final ArrayList<T> array;

		private ArrayData(ArrayList<T> array) {
			this.array = array;
		}
		@Override public T get(int i) {
			return array.get(i);
		}
	}
	/** One link in a diff chain. exactly the same as 'oldData', except the 'diffIdx'th entry is 'diffVal' */
	private static class DiffData<T> implements Data<T> {
		private final int diffIdx;
		private T diffVal;
		private PersistentArrayList<T> oldData;

		private DiffData(int diffIdx, T diffVal, PersistentArrayList<T> oldData) {
			this.diffIdx = diffIdx;
			this.diffVal = diffVal;
			this.oldData = oldData;
		}
		@Override public T get(int i) {
			return (i == diffIdx) ? diffVal : oldData.get(i);
		}
	}

	private PersistentArrayList(Data<T> data, int size) {
		this.data = data;
		this.size = size;
	}

	/** Creates a new PersistentArrayList with a copy of the given collection */
	public PersistentArrayList(Collection<T> collection) {
		final ArrayList<T> list = Lists.newArrayList(collection); // defensive copy
		this.data = new ArrayData<T>(list);
		this.size = list.size();
	}

	/** Gets the value stored at index i. Note that this triggers a reroot. */
	public T get(int i) {
		checkElementIndex(i, size);
		return data.get(i);
	}

	/**
	 * Gets a PersistentArrayList which is the same as this, except it has 'val' at position 'i'.
	 * The returned PersistentArrayList will be backed by an array, and this will be a diff pointing to it.
	 * @param i the index to update
	 * @param val the value to put at index 'i'
	 * @return a new PersistentArrayList which is the same as this, except it has 'val' at position 'i'.
	 */
	public PersistentArrayList<T> set(int i, T val) {
		checkElementIndex(i, size);
		if (get(i) == val) return this;
		reRoot();
		final ArrayData<T> arrayData = (ArrayData<T>) data;
		final T oldVal = arrayData.array.get(i);
		arrayData.array.set(i, val);
		final PersistentArrayList<T> newArrayList = new PersistentArrayList<T>(data, size);
		data = new DiffData<T>(i, oldVal, newArrayList);
		return newArrayList;
	}

	@Override public int size() {
		return size;
	}

	/**
	 * Reverses all of the diffs on our chain so that this.data becomes an ArrayData.
	 * We do this to make the most recently accessed PersistentArrayList the most efficient,
	 * under the assumption that it is more likely to be accessed again in the future.
	 * Afterward, this.data is guaranteed to be of type ArrayData.
	 */
	private void reRoot() {
		if (data instanceof DiffData) {
			// recursively reroot oldData
			final DiffData<T> diffData = (DiffData<T>) data;
			diffData.oldData.reRoot();
			// swap oldData with this, reversing the diff
			final ArrayData<T> rootedArray = (ArrayData<T>) diffData.oldData.data;
			final T oldVal = rootedArray.array.get(diffData.diffIdx);
			rootedArray.array.set(diffData.diffIdx, diffData.diffVal);
			diffData.oldData.data = new DiffData<T>(diffData.diffIdx, oldVal, this);
			data = rootedArray;
		}
	}

	@Override public Iterator<T> iterator() {
		reRoot();
		return ((ArrayData<T>) data).array.iterator();
	}
}

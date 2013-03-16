package edu.cmu.cs.lti.ark.cle;

import com.google.common.collect.Lists;
import org.junit.Test;

import static edu.cmu.cs.lti.ark.cle.PersistentArrayList.copyOf;
import static org.junit.Assert.assertEquals;

/**
 * @author sthomson@cs.cmu.edu
 */
public class PersistentArrayListTest {
	@Test
	public void testPersistence() {
		final PersistentArrayList<String> strings1 = copyOf(Lists.newArrayList("a", "b", "c"));
		final PersistentArrayList<String> strings2 = strings1.set(1, "b2");
		final PersistentArrayList<String> strings3 = strings1.set(2, "c3");
		final PersistentArrayList<String> strings4 = strings3.set(0, "a4");
		assertEquals("a", strings1.get(0));
		assertEquals("a", strings2.get(0));
		assertEquals("a", strings3.get(0));
		assertEquals("a4", strings4.get(0));
		assertEquals("b", strings1.get(1));
		assertEquals("b2", strings2.get(1));
		assertEquals("b", strings3.get(1));
		assertEquals("b", strings4.get(1));
		assertEquals("c", strings1.get(2));
		assertEquals("c", strings2.get(2));
		assertEquals("c3", strings3.get(2));
		assertEquals("c3", strings4.get(2));
	}

}
